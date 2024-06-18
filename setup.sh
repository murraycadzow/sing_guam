# $1 is the account list file

for account in $(cat $1)
do
	echo ${account}
done
