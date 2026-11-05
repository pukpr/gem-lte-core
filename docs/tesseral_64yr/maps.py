"""Digitize the coloured gauge dots of Fig. 3 (a) amplitude and (d) phase maps into (lat, lon, amp, phase)."""
import sys
from PIL import Image; import numpy as np, json
from scipy import ndimage
im=np.array(Image.open(sys.argv[1] if len(sys.argv) > 1 else 'fig3.png').convert('RGB')).astype(float)
sat=(im.max(2)-im.min(2))
lonf=lambda x: ((x-68.0)/84.3*90.0)
latf=lambda y: (129.5-y)/56.25*60.0
cba=im[225,118:357]; vxa=(np.arange(118,357)-117)/46.6*6
cbd=im[448,118:357]; vxd=(np.arange(118,357)-237)/59.75*90
def val(rgb,cb,vx):
    d=((cb-rgb)**2).sum(1); i=d.argmin(); return vx[i], np.sqrt(d[i])
m=(sat>70); m[:46,:]=False; m[215:,:]=False; m[:,:64]=False; m[:,418:]=False
lab,n=ndimage.label(m,structure=np.ones((3,3)))
pts=[]
for i in range(1,n+1):
    ys,xs=np.where(lab==i)
    if len(xs)<6: continue
    # split big blobs into ~6x6 cells so clustered gauges keep their own colours
    for cy in range(ys.min(),ys.max()+1,6):
        for cx in range(xs.min(),xs.max()+1,6):
            s=(ys>=cy)&(ys<cy+6)&(xs>=cx)&(xs<cx+6)
            if s.sum()<6: continue
            yy,xx=ys[s],xs[s]
            ca=im[yy,xx].mean(0); cd=im[yy+223,xx].mean(0)
            a,ea=val(ca,cba,vxa); p,ep=val(cd,cbd,vxd)
            if sat[yy+223,xx].mean()<50: continue
            pts.append(dict(lat=float(latf(yy.mean())),lon=float(lonf(xx.mean())),amp=float(a),phase=float(p),err_a=float(ea),err_p=float(ep)))
json.dump(pts,open('map_points.json','w'),indent=0)
print(len(pts),'cells'); 
import collections
for q in pts[:5]: print(q)
