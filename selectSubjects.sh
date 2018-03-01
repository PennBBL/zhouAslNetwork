for i in `cat commonSubjects.txt`
do
grep "$i" cbfNetworkList.txt
done