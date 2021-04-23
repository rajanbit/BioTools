cd ..
python scripts/gc_analyzer.py -n bin/$1 -o temp/xox333.tsv -c 75000  
python scripts/gc_analyzer.py -n bin/$2 -o temp/x1x333.tsv -c 75000 

if [[ "stat -c%s bin/$1" > "stat -c%s bin/$2" ]];then
paste temp/x1x333.tsv temp/xox333.tsv | awk '{print $1 "\t" $2 "\t" $4 }' > temp/xox88.tsv
else
paste temp/xox333.tsv temp/x1x333.tsv | awk '{print $1 "\t" $2 "\t" $4 }' > temp/xox88.tsv
fi

sed -i '1s/^/Pos \t GC1 \t GC2 \n/' temp/xox88.tsv
Rscript scripts/plotter.R temp/xox88.tsv

