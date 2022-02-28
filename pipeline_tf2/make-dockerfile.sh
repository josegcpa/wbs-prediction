cat Dockerfile-header > Dockerfile

x=''
i=0
M=1000
for p in $(cat requirements.txt | grep -v "[jJ]upyter" | grep -v "jedi" | grep -v "ipyt" | grep -v "Mark" | grep -v "Jinja" | grep -v "ipyker" | grep -v "^nb" | grep -v "noteb" | grep -v "tensorboa")
#for p in $(cat requirements.txt)
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

if [[ ${#x} -gt 0 ]]
then
    echo RUN pip3 install $x >> Dockerfile
fi

echo "" >> Dockerfile

echo "RUN apt-get update -y && apt-get install -y libgl1-mesa-dev" >> Dockerfile
echo "" >> Dockerfile
echo COPY scripts scripts >> Dockerfile
echo COPY Snakefile Snakefile >> Dockerfile
echo COPY pipeline.sh run-slide >> Dockerfile
echo COPY run-folder.sh run-folder >> Dockerfile
echo COPY test-packages.py test-packages.py >> Dockerfile
echo "" >> Dockerfile
echo RUN chmod +x /app/run-slide >> Dockerfile
echo RUN chmod +x /app/run-folder >> Dockerfile
