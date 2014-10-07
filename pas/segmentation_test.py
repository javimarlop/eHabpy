
i.segment gr=segm thr=0.5 out=sn_segm_vars052 meth='region_growing' simil='euclidean' minsi=100 iter=20 --o
r.pintame sn_segm_vars052
r.mask vect=sn
r.mapcalc 'sn_segm_vars052m = sn_segm_vars052' --o
g.remove MASK
r.pintame sn_segm_vars052m
i.group gr=segm2 in=sn_segm_vars052m
i.segment gr=segm2 thr=0.01 out=sn_segm_vars052pa meth='region_growing' simil='euclidean' minsi=100 iter=20 --o
r.pintame sn_segm_vars052pa
r.stats -nc sn_segm_vars052pa sort=desc
r.mapcalc 'sn_segm_vars052paf = nmode(sn_segm_vars052pa[-1,0],sn_segm_vars052pa[-1,1],sn_segm_vars052pa[-1,-1],sn_segm_vars052pa[0,1],sn_segm_vars052pa[0,-1],sn_segm_vars052pa[1,0],sn_segm_vars052pa[1,1],sn_segm_vars052pa[1,-1])' --o
r.pintame sn_segm_vars052paf
r.stats -nc sn_segm_vars052paf sort=desc
r.to.vect -v input=sn_segm_vars052paf out=sn_segm_vars052paf type='area' --o
#v.clean input=sn_segm_vars052pa out=sn_segm_vars052paf tool='rmarea' thres=100000 --o
r.neighbors in=sn_segm_vars052pa out=sn_segm_vars052paf_ecotones meth=div size=3 --o # apply this to obtain areas acting as ecotones!
#d.vect -c sn_segm_vars052paf


