for f in proc*.txt;
do sort $f|uniq; 
done|sort|uniq -c -d|sort -nr -k1,1 > commonSubjects.txt