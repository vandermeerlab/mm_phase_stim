import numpy as np
import cv2
from sklearn.decomposition import PCA
from scipy.stats import linregress
from scipy.signal import find_peaks

# very important to crop so that can only see fur, no feet or other objects

lk_params = dict( winSize  = (15, 15),
                  maxLevel = 2,
                  criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 0.03))

feature_params = dict( maxCorners = 500,
                       qualityLevel = 0.3,
                       minDistance = 7,
                       blockSize = 7 )


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return np.isnan(y), lambda z: z.nonzero()[0]


class VideographicPlethysmography:
	def __init__(self, video_src,fps,targetNumPts=100,startX=0,endX=None,startY=0,endY=None,view=False):
		self.video_src = video_src
		self.targetNumPts = targetNumPts
		self.startX = startX
		self.endX = endX
		self.startY = startY
		self.endY = endY
		self.view=view
		self.ptList = [] # define lists that will hold point data as [frameNum, id, x position, yposition]
		self.currPtIDs = [] # define list that will hold the current pt IDs so we can keep track of point identity
		self.lastAddedPtID = -1; # remembers how many points we have added so we can add a unique ID each time a new point is found
		self.cam = cv2.VideoCapture(video_src) # initialize capture device, in this case a file
		self.frameNum = 0
		self.fps = fps
		self.pca = PCA(n_components=1)
		self.motion_energy = []
	def track_vid(self):
		_ret, frame = self.cam.read() # read our first frame
		if type(self.endX) == type(None) or self.endX > frame.shape[0]:
			self.endX = frame.shape[0]
		if type(self.endY) == type(None) or self.endY > frame.shape[1]:
			self.endY = frame.shape[1]
		frame_gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)[self.startX:self.endX,self.startY:self.endY] # convert to grayscale
		mask = np.zeros_like(frame_gray) # create mask for finding points
		mask[:] = 255
		newPtCoords = cv2.goodFeaturesToTrack(frame_gray, mask = mask, **feature_params) # find good points to track
		prev_gray = frame_gray # set as previous frame to set up optic flow calculation
		# add newly found pts to point list, update point id info accordingly
		for newPt in range(len(newPtCoords)):
			self.ptList.append([self.frameNum, self.lastAddedPtID+1,newPtCoords[newPt,0,0],newPtCoords[newPt,0,1]]) # [frame, id, x position, yposition]
			self.currPtIDs.append(self.lastAddedPtID+1)
			self.lastAddedPtID += 1
		p = newPtCoords; # set up points for propagation
		# continuously propagate points from previous step and find new points until no more frames found
		while self.cam.isOpened():
			ret, frame = self.cam.read()
			if ret:
				self.frameNum += 1
				frame_gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)[self.startX:self.endX,self.startY:self.endY]
				self.motion_energy.append(np.sum(np.abs(frame_gray-prev_gray)))
				img0, img1 = prev_gray, frame_gray
				p0 = p
				p1, _st, _err = cv2.calcOpticalFlowPyrLK(img0, img1, p0.astype('float32'), None, **lk_params) # propagate forward
				if type(p1) != type(None): # if backpropagation doesnt catastrophically fails 
					p0r, _st, _err = cv2.calcOpticalFlowPyrLK(img1, img0, p1.astype('float32'), None, **lk_params) # propagate backward
					d = abs(p0-p0r).reshape(-1, 2).max(-1) # take error in backward and forward error
					good = d < 1 # threshold error
				else: # catastrophic failure, stop tracking previous points
					good = p0>np.inf
				keepPTIDs = []
				for i in range(len(good)): # make new list for ids that were successfully propagated and append the propagated points to list
					if good[i]:
						keepPTIDs.append(self.currPtIDs[i])
						self.ptList.append([self.frameNum,self.currPtIDs[i],p1[i,0,0],p1[i,0,1]]) # add successfully propagated points to ptList
				self.currPtIDs=keepPTIDs
				newPtCoords = None
				if len(self.currPtIDs) < self.targetNumPts:
					if type(p1) != type(None):
						blockPoints =  np.squeeze(np.round(p1).astype(int)) # format found points
					else:
						blockPoints = np.zeros(1) # give it shape that will prevent any masking
					mask = np.zeros_like(frame_gray)
					mask[:] = 255
					if len(blockPoints.shape)>1:
						blockPoints[blockPoints[:,1]>=frame_gray.shape[0],1] = frame_gray.shape[0]-1
						blockPoints[blockPoints[:,1]<0] = 0
						blockPoints[blockPoints[:,0]>=frame_gray.shape[1],0] = frame_gray.shape[1]-1
						blockPoints[blockPoints[:,0]<0] = 0
						mask[blockPoints[:,1],blockPoints[:,0]] = 0 # block out already found points, prevents duplicates
					newPtCoords = cv2.goodFeaturesToTrack(frame_gray, mask = mask, **feature_params) # get newly found points
				numNewPts = 0
				if type(newPtCoords) != type(None):
					for newPt in range(len(newPtCoords)): # for each new point, add it to pt list and update point id info 
						self.ptList.append([self.frameNum, self.lastAddedPtID+1,newPtCoords[newPt,0,0],newPtCoords[newPt,0,1]]) # [frame, id, x position, yposition]
						self.currPtIDs.append(self.lastAddedPtID+1)
						self.lastAddedPtID += 1
					numNewPts,_,_ = np.shape(newPtCoords) # get number of newly found points
				numOldPts = np.sum(good) # get number of kept old points
				p = np.zeros((numOldPts+numNewPts,1,2)) # create array for new point data to be propagated in next step
				added = 0
				for i in range(len(good)): # add old points to array
					if good[i]:
						p[added,0,:] = p1[i,0,:] # added was just i, changed to account for some old points that are no longer tracked
						added+=1
				if type(newPtCoords) != (None):
					p[numOldPts:,:,:] = newPtCoords # add new points to array
				prev_gray = frame_gray # set up optic flow algorithm
				if self.view:
					real_vis = frame.copy()[self.startX:self.endX,self.startY:self.endY]
					for i in range(p.shape[0]):
						circ = cv2.circle(real_vis, (int(p[i,0,0]), int(p[i,0,1])), 2, (0, 255, 0), -1)
					cv2.imshow('lk_track',real_vis)
					if cv2.waitKey(1) & 0xFF == ord('q'):
						break
			else:
				self.cam.release()
				break
	def plethysmography(self):
		ptsByID = [] 
		for i in range(self.lastAddedPtID+1): # add list for every point
			ptsByID.append([])
		for i in self.ptList:
			ptsByID[i[1]].append([i[0],i[2],i[3]]) # add frame and x y coord in id specific list
		pcaList = []
		for i in ptsByID:
			if len(i) > self.fps: # only compute for points we could track for a second
				data = np.zeros((len(i),2))
				for ii in range(len(i)):
					data[ii,0] = i[ii][1]
					data[ii,1] = i[ii][2]
				transform = self.pca.fit_transform(data) # project movement onto principal axis of motion
				linearModel = linregress(np.cumsum(np.ones(np.size(transform))),np.squeeze(transform))
				detrended = np.squeeze(transform)-(np.cumsum(np.ones(np.size(transform)))*linearModel.slope+linearModel.intercept) # account for drift 
				pcaList.append([i[0][0],i[-1][0],detrended]) # append to pcaList
			else:
				pcaList.append([-1,-1,-1])
		self.meanPCAProj = np.zeros((self.frameNum))
		self.stdPCAProj = np.copy(self.meanPCAProj)
		self.nPCAProj = np.copy(self.meanPCAProj)
		for i in range(self.frameNum):
			pcaVals = []
			for ii in pcaList:
				if i >= ii[0] and i <= ii[1]:
					pcaVals.append(ii[2][i-ii[0]])
			quickData = np.array(pcaVals)
			try: # modify to check if empty array
				self.meanPCAProj[i] = np.mean(quickData)
				self.stdPCAProj[i] = np.std(quickData)
				self.nPCAProj[i] = len(pcaVals)
			except:
				self.meanPCAProj[i] = np.nan
				self.stdPCAProj[i] = np.nan
				self.nPCAProj[i] = np.nan
		nans, x= nan_helper(self.meanPCAProj)
		self.meanPCAProj[nans]= np.interp(x(nans), x(~nans), self.meanPCAProj[~nans])
		sP = moving_average(self.meanPCAProj,5)
		inhalePks,_ = find_peaks(sP,prominence=(0.5,None))
		exhalePks,_ = find_peaks(-sP,prominence=(0.5,None))
		numInhalePks = inhalePks.shape[0]
		numExhalePks = exhalePks.shape[0]
		exhaleStart = exhalePks[0]<inhalePks[0]
		self.instantRespiratoryRate = np.zeros(np.shape(self.meanPCAProj))*np.nan
		self.instantRespiratoryAmp = np.zeros(np.shape(self.meanPCAProj))*np.nan
		for i in range(numInhalePks-1):
			self.instantRespiratoryRate[inhalePks[i]:inhalePks[i+1]] =  1/((inhalePks[i+1]-inhalePks[i])/self.fps)
		if exhaleStart:
			for i in range(min(numInhalePks,numExhalePks)-1):
				self.instantRespiratoryAmp[inhalePks[i]:inhalePks[i+1]] = sP[inhalePks[i]]-sP[exhalePks[i]]
		else:
			for i in range(1,min(numInhalePks,numExhalePks)-1):
				self.instantRespiratoryAmp[inhalePks[i]:inhalePks[i+1]] = sP[inhalePks[i]]-sP[exhalePks[i-1]]
	def getPCByPoints(self):
		ptsByID = [] 
		for i in range(self.lastAddedPtID+1): # add list for every point
			ptsByID.append([])
		for i in self.ptList:
			ptsByID[i[1]].append([i[0],i[2],i[3]]) # add frame and x y coord in id specific list
		pcaList = []
		for i in ptsByID:
			if len(i) > self.fps: # only compute for points we could track for a second
				data = np.zeros((len(i),2))
				for ii in range(len(i)):
					data[ii,0] = i[ii][1]
					data[ii,1] = i[ii][2]
				transform = self.pca.fit_transform(data) # project movement onto principal axis of motion
				linearModel = linregress(np.cumsum(np.ones(np.size(transform))),np.squeeze(transform))
				detrended = np.squeeze(transform)-(np.cumsum(np.ones(np.size(transform)))*linearModel.slope+linearModel.intercept) # account for drift 
				pcaList.append([i[0][0],i[-1][0],detrended]) # append to pcaList
			else:
				pcaList.append([-1,-1,-1])
		dataMat = np.zeros((len(pcaList),self.frameNum+1))
		for i in range(len(pcaList)):
			dataMat[i,pcaList[i][0]:pcaList[i][1]+1] = pcaList[i][2]
		return dataMat
	def getInstantRespiratoryRate(self):
		return self.instantRespiratoryRate
	def getInstantRespiratoryAmp(self):
		return self.instantRespiratoryAmp
	def getMeanPCAProj(self):
		return self.meanPCAProj
	def getStdPCAProj(self):
		return self.getStdPCAProj
	def getNumTrackedPoints(self):
		return self.nPCAProj
	def getFirstFrame(self):
		self.cam = cv2.VideoCapture(self.video_src)
		ret, frame = self.cam.read()
		return frame.copy()
	def getMotionEnergy(self):
		return self.motion_energy





