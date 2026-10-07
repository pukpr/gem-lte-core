#!/usr/bin/env python3
"""Compare the ZONAL fit (user's runs/f_amo_ZONAL, snapshot) with the TABLE fit (savebug/w4, CC 0.833):
model output, composite winding sum on each fit's own manifold, and the manifold itself.
For each: CC at lag 0, best lag (+-60 months), and the residual after the best lag and an affine
rescale (its variance share and the periods carrying it)."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z=np.loadtxt("zonal_lte_results.csv",delimiter=","); T=np.loadtxt("table_lte_results.csv",delimiter=",")
n=min(len(Z),len(T)); Z,T=Z[:n],T[:n]; assert np.allclose(Z[:,0],T[:,0]); t=Z[:,0]; data=Z[:,2]
def comp(res, w):   # composite winding sum on the fit's own manifold (column 4)
    F=res[:,3]; j=json.load(open(w))["k_amp_phase"]
    return sum(a*np.sin(2*np.pi*k*F+ph) for k,a,ph in j if a!=0)
S={"model output (col 2)":(Z[:,1],T[:,1]),"composite winding sum":(comp(Z,"zonal_windings.json"),comp(T,"table_windings.json")),"manifold F2 (col 4)":(Z[:,3],T[:,3])}
def lagcc(a,b,L):   # CC of a(t) with b(t-L)
    if L>=0: return np.corrcoef(a[L:],b[:n-L])[0,1]
    return np.corrcoef(a[:n+L],b[-L:])[0,1]
f=np.fft.rfftfreq(n,1/12)
def bands(r):
    f=np.fft.rfftfreq(len(r),1/12); P=np.abs(np.fft.rfft(r-r.mean()))**2; tot=P[1:].sum()
    return {lab:100*P[(f>=lo)&(f<hi)].sum()/tot for lab,lo,hi in (("<1 yr",1,6.1),("1-2 yr",0.5,1),("2-7 yr",1/7,0.5),("7-20 yr",0.05,1/7),(">20 yr",0.0001,0.05))}
print(f"data CC: ZONAL model {np.corrcoef(Z[:,1],data)[0,1]:.4f}, TABLE model {np.corrcoef(T[:,1],data)[0,1]:.4f}")
fig,ax=plt.subplots(3,2,figsize=(14,10),gridspec_kw=dict(width_ratios=[3,1],hspace=0.35))
for row,(lab,(a,b)) in enumerate(S.items()):
    lags=np.arange(-60,61); c=np.array([lagcc(a,b,L) for L in lags]); L=lags[np.argmax(np.abs(c))]
    bs=b[:n-L] if L>=0 else b[-L:]; as_=a[L:] if L>=0 else a[:n+L]; ts=t[L:] if L>=0 else t[:n+L]
    X=np.column_stack([bs,np.ones_like(bs)]); cf,*_=np.linalg.lstsq(X,as_,rcond=None); res=as_-X@cf
    bd=bands(res)
    print(f"{lab}: CC lag0 {c[60]:+.4f}; best lag {L:+d} months (ZONAL lags TABLE by {L} mo) CC {c[60+L]:+.4f}; residual after shift+scale = {100*res.var()/as_.var():.1f}% of ZONAL variance; residual bands: " + ", ".join(f"{k} {v:.0f}%" for k,v in bd.items()))
    fr=np.fft.rfftfreq(len(res),1/12); P=np.abs(np.fft.rfft(res-res.mean()))**2; top=np.argsort(P[1:])[-4:][::-1]+1
    print(f"   residual's strongest periods: " + ", ".join(f"{1/fr[k]:.2f} yr" for k in top))
    ax[row,0].plot(t,(a-a.mean())/a.std(),lw=0.9,color="#a05fd0",label="ZONAL"); ax[row,0].plot(t,(b-b.mean())/b.std(),lw=0.9,color="#eb6834",label="TABLE",alpha=0.8)
    ax[row,0].plot(ts,(res-res.mean())/as_.std()-3.5,lw=0.7,color="#555",label=f"residual (after lag {L:+d} mo + scale), offset")
    ax[row,0].set_title(f"{lab}: CC lag0 {c[60]:+.3f}, best lag {L:+d} mo CC {c[60+L]:+.3f}",loc="left",fontsize=10); ax[row,0].legend(frameon=False,fontsize=7,ncol=3,loc="upper left"); ax[row,0].grid(alpha=.3)
    ax[row,1].plot(lags,c,color="#2a78d6"); ax[row,1].axvline(0,color="#c3c2b7"); ax[row,1].set_title("cross-correlation vs lag (months)",fontsize=9); ax[row,1].grid(alpha=.3)
ax[2,0].set_xlabel("Year"); fig.suptitle("AMO: ZONAL (astronomical forcing) vs TABLE (constituent table) fits",x=0.06,ha="left")
fig.savefig("zonal_vs_table.png",dpi=110,bbox_inches="tight"); print("plot: compare/zonal_vs_table.png")
