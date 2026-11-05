"""Digitize the black estimation dots of Fig. 3 (b, c, e, f) of the 64-yr sea-level paper."""
import sys
from PIL import Image; import numpy as np, json
from scipy import ndimage
im=np.array(Image.open(sys.argv[1] if len(sys.argv) > 1 else 'fig3.png').convert('RGB')).astype(int)
blk=(im.max(2)<110)&((im[:,:,2]-im[:,:,0])<40)
def lin(px, val):
    a,b=np.polyfit(px,val,1); return lambda p: a*p+b
cal={'b':(lin([449,516,583,650.5],[-30,0,30,60]), lin([70,118.5,166,214.5],[24,16,8,0])),
     'c':(lin([744.5,814.5,882,950.5],[-90,0,90,180]), lin([70,118.5,166,214.5],[24,16,8,0])),
     'e':(lin([449,516,583,650.5],[-30,0,30,60]), lin([268.5,311,353,395.5,437.5],[180,90,0,-90,-180])),
     'f':(lin([744.5,814.5,883,951],[-90,0,90,180]), lin([268.5,311,353,395.5,437.5],[180,90,0,-90,-180]))}
frames={'b':(420,673,46,214),'c':(716,969,46,214),'e':(420,673,268,437),'f':(716,969,268,437)}
excl={'b':[(420,460,46,72),(588,673,184,214)],'c':[(716,755,46,72)],'e':[(420,460,268,294)],'f':[(716,755,268,294),(898,969,404,437)]}
out={}
for k,(x0,x1,y0,y1) in frames.items():
    m=np.zeros_like(blk); m[y0+2:y1-1,x0+2:x1-1]=blk[y0+2:y1-1,x0+2:x1-1]
    for (a,b,c,d) in excl[k]: m[c:d,a:b]=False
    lab,n=ndimage.label(m,structure=np.ones((3,3)))
    pts=[]
    for i in range(1,n+1):
        ys,xs=np.where(lab==i); area=len(xs)
        if area<4: continue
        cnt=max(1,int(round(area/14)))   # ~14 px per dot; merged blobs count multiple
        fx,fy=cal[k]
        pts+= [(float(fx(xs.mean())),float(fy(ys.mean())))]*cnt
    out[k]=pts; print(k,len(pts))
json.dump(out,open('fig3_points.json','w'))
