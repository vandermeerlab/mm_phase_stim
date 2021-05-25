import VideographicPlethysmography as vp
import matplotlib.pyplot as plt
vpObj = vp.VideographicPlethysmography('./IMG_3569.mov',fps=60,startX=300,endX=-150,startY=250,endY=-300,view=True)
vpObj.track_vid()
dataMat = vpObj.getPCByPoints()
plt.matshow(dataMat)