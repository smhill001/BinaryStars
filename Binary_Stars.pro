x=read_bmp('C:\Astronomy\Images\2010\01\30\EtaCas_5500mm_120sec__Rot1_0000_DarkStack600Ref50Gamma2x.bmp')
y=reform(x(1,*,*),1280,960)
xtv,y
cntrd,y,630,495,xcen,ycen
5
print,xcen,ycen