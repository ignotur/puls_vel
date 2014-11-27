#!/bin/bash



scp -i ~/.ssh/id_rsa -r aigoshev@coma.science.ru.nl:/home/aigoshev/puls_vel/analysis/profiles .
cp profiles/velocity/*.dat profiles/

## Copy programs which find the best parameter
cp ../../brisken/analyse_model_1D_maxw.out .
cp ../../brisken/analyse_model_1D_norm.out .
cp ../../brisken/analyse_model_1D_pach.out .
cp ../../brisken/analyse_model_1D_uni.out .
cp ../../brisken/analyse_model_1D_sum_maxw.out .

## Copy program which estimets credential interval
cp ../../brisken/credential_1D_maxw.out .
cp ../../brisken/credential_1D_norm.out .
cp ../../brisken/credential_1D_pach.out .
cp ../../brisken/credential_1D_sum_maxw.out .

## Analyse the Maxwellian distribution
./analyse_model_1D_maxw.out
tail -n 2 model.dat > tmp
param=`head -n 1 tmp`
./credential_1D_maxw.out $param

## Analyse the normal distribution
./analyse_model_1D_norm.out
tail -n 2 model.dat > tmp
param=`head -n 1 tmp`
./credential_1D_norm.out $param

## Analyse the Paczynski distribution
./analyse_model_1D_pach.out
tail -n 2 model.dat > tmp
param=`head -n 1 tmp`
./credential_1D_pach.out $param

## Analyse the Uniform distribution
./analyse_model_1D_pach.out
echo "*" >> model.dat
echo "*" >> model.dat

## Analyse the sum of two Maxwellian distribution
./analyse_model_1D_sum_maxw.out
tail -n 2 model.dat > tmp
param=`head -n 1 tmp`
./credential_1D_sum_maxw.out $param

## Non-parametric study
cp ../standard_28092014/non-parametric/TC93thick/1D_joint/fy/skr_obr .
cp ../standard_28092014/non-parametric/TC93thick/1D_joint/fy/sum_profile.out .
cp ../standard_28092014/non-parametric/TC93thick/1D_joint/fy/monoton.py .
cp ../standard_28092014/non-parametric/TC93thick/1D_joint/fy/1D_to_cumul_distr.out .

./skr_obr



## make a .pdf file with velocity and distance PDFs
## First distances
cd profiles/distance/
cp ../../../../brisken/develop/tests/newer_test_26112014/plot_dist1.gnu  .
#scp aigoshev@coma.science.ru.nl:/home/aigoshev/puls_vel/catalogue/what_is_name.py .
cp ../../../what_is_name.py .
scp -i ~/.ssh/id_rsa  aigoshev@coma.science.ru.nl:/home/aigoshev/puls_vel/catalogue/name_list .

num=`ls -l profile_*.dat | wc -l`


for (( i = 0; i < $num; i++))
do
cp profile_${i}.dat tmp
name=`python what_is_name.py $i`
gnuplot -e "filename='${name}'" plot_dist1.gnu
convert -density 150 tmp.eps profile_${i}.pdf
done
pdftk profile_*.pdf cat output distance_pdf.pdf

## Second velocity
cd ../../

cd profiles/velocity

cp ../../../standard_28092014/parametric/TC93thick/1D_joint/fy/plot_overall_vel.gnu .
cp ../../profiles/distance/what_is_name.py .
cp ../../profiles/distance/name_list .

num=`ls -l profile_*.dat | wc -l`

for (( i = 0; i < $num; i++ ))
do
name=`python what_is_name.py $i`
gnuplot -e "filename='../../../analysis_done_20112014/profiles/velocity/profile_${i}.dat'; filename2='profile_${i}.dat'; filename3='${name}'" plot_overall_vel.gnu
convert -density 150 tmp.eps profile_${i}.pdf
done
pdftk profile_*.pdf cat output velocity_pdf.pdf

