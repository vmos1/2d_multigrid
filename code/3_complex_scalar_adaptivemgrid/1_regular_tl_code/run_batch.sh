L=64
blk=2
ndof=2
for lev in 4 3 2 1 0;
do for m in 0.0 -0.02 -0.03 -0.04;
    do echo $L,3, $blk, 0, $m,$lev;
    ./a1 $L 3 $blk 0 $m $lev;
#    echo results_phi.txt results_phi_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_ndof"$ndof".txt
    mv results_phi.txt results_phi_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_ndof"$ndof".txt
    mv results_residue.txt results_residue_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_ndof"$ndof".txt
    done;
done
