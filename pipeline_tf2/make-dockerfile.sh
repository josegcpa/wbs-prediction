cat Dockerfile-header > Dockerfile

x=''
i=0
M=10
for p in $(cat requirements.txt)
do
    x=$(echo $x $p)
    i=$(($i+1))
    if [[ $i == $M ]]
    then
        echo RUN pip3 install $x >> Dockerfile
        x=''
        i=0
    fi
done
echo RUN pip3 install $x >> Dockerfile
