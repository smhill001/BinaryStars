pro doublestar,tiffname

;doublestar,'C:\Astronomy\Images\2008\04\14\XIUMALVA.TIF'
;xrng=[450,500]
;yrng=[225,275]

;doublestar,'C:\Astronomy\Images\2008\02\19\CASTRLVA.TIF'
yrng=[0,1]
xrng=[0,1]

x=read_tiff(tiffname)
y=reform(x(1,*,*),1800,1200)

;xtv,y
;stop

y=float(y)
z1=total(y(*,xrng(0):xrng(1)),2)
print,where(z1 eq max(z1(yrng(0):yrng(1))))
;stop

z2=total(y(yrng(0):yrng(1),*),1)
print,where(z2 eq max(z2(xrng(0):xrng(1))))

stop
end

