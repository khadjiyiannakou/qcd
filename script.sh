paste three_up_axial_1.dat three_down_axial_1.dat > temp1.dat
paste three_up_axial_2.dat three_down_axial_2.dat > temp2.dat
paste three_up_axial_3.dat three_down_axial_3.dat > temp3.dat

awk <temp1.dat '{printf("%d %d %e %e\n",$1,$2,($3-$7)*0.5,($4-$7)*0.5) }' > kale1.dat
awk <temp2.dat '{printf("%d %d %e %e\n",$1,$2,($3-$7)*0.5,($4-$7)*0.5) }' > kale2.dat
awk <temp3.dat '{printf("%d %d %e %e\n",$1,$2,($3-$7)*0.5,($4-$7)*0.5) }' > kale3.dat

paste temp1.dat temp2.dat temp3.dat > temp.dat

awk <temp.dat '{printf("%d %d %e %e\n",$1,$2,$3+$7+$11,$4+$8+$12) }' > final.dat

sed -n "1,11p" final.dat > final_zero.dat

awk < final_zero.dat '{printf("%d %d %f\n",$1,$2,$3/(1.045519e-10))}' > ratio.dat

