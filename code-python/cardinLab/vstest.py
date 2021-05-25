import VideographicPlethysmography as vp
import matplotlib.pyplot as plt
vpObj = vp.VideographicPlethysmography('./test2.mp4',fps=50,startX=123,endX=250,startY=74,endY=430,view=True)
vpObj.track_vid()
dataMat = vpObj.getPCByPoints()
plt.matshow(dataMat)
plt.show()