#!/bin/bash
#db_dir=/Users/mapostolides/Drive/fusion-pipeline/trusight
db_dir=/Users/mapostolides/MetaFusion.clinical/reference_files
if [ "$#" -ne 1 ]; then
    echo "You must choose either "add" or "remove" "
	echo USAGE: sh update_false_positives.sh add
	echo USAGE: sh update_false_positives.sh remove
	echo USAGE: sh update_false_positives.sh view
	exit 1
fi

# SELECT MODE: ADD OR REMOVE
if [ $1 == add ];then 
	mode=add;
    echo ADDING FUSION TO DATABASE
elif [ $1 == remove ]; then 
	mode=remove; 
    echo REMOVING FUSION FROM DATABASE 
elif [ $1 == view ]; then
	mode=view;
	echo VIEWING DATABASE
else
	echo USAGE: sh update_false_positives.sh add
	echo USAGE: sh update_false_positives.sh remove
	echo USAGE: sh update_false_positives.sh view
	exit 1
fi

if [ $mode == add ]; then

echo gene1:
read gene1
echo gene2:
read gene2
    result=$(sqlite3 $db_dir/historical_fusions.db "SELECT * FROM false_positives WHERE gene1='$gene1' AND gene2='$gene2';")
	if [ -z "$result" ]
	then
    	  echo FUSION \"$gene1"--"$gene2\" NOT YET IN DATABASE. ADDING 
		  sqlite3 $db_dir/historical_fusions.db "INSERT INTO false_positives VALUES('$gene1','$gene2');select * from false_positives;"
	else
    	 echo $result ALREADY IN DATABASE 
	fi
	#sqlite3 $db_dir/historical_fusions.db "INSERT INTO false_positives VALUES('$gene1','$gene2');select * from false_positives;"

elif [ $mode == remove ]; then
echo gene1:
read gene1
echo gene2:
read gene2
    result=$(sqlite3 $db_dir/historical_fusions.db "SELECT * FROM false_positives WHERE gene1='$gene1' AND gene2='$gene2';")
	    if [ -z "$result" ]
    then
          echo FUSION \"$gene1"--"$gene2\" NOT IN DATABASE. CANNOT REMOVE 
    else
         echo $result IS IN DATABASE. REMOVING
		 sqlite3 $db_dir/historical_fusions.db "DELETE FROM false_positives WHERE gene1='$gene1' AND gene2='$gene2';select * from false_positives;"
    fi
else
	sqlite3 $db_dir/historical_fusions.db "select * from false_positives;"
fi
