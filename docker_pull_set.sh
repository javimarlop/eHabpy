### Javier Martinez-Lopez
# UTF-8
# 2014

docker pull javimarlop/ubuntugis-docker # clone docker image
#docker run --rm -tiv /home/webuser/grassdb:/opt/grassdb javimarlop/ubuntugis-docker # with shared folder
docker run --rm --privileged -ti javimarlop/ubuntugis-docker /bin/bash

useradd -m --uid xxx majavie
usermod -a -G fuse majavie
usermod -s /bin/bash majavie

su majavie
cd /home/majavie
mkdir /home/majavie/hanksgrass7
mkdir /home/majavie/data
sshfs xxx@yyy:/local1/majavie/hanksgrass7 /home/majavie/hanksgrass7/ # mount grass db
sshfs xxx@yyy:/local0/majavie/eHabpy /home/majavie/data/ # mount ehabpy dir with all data
#fusermount -u /home/webuser/hanksgrass7/
#cd data/pas
python docker_test_grass.py

## tzar localrun
#wget https://tzar-framework.atlassian.net/wiki/download/attachments/4980739/tzar-0.5.4.jar?version=1&modificationDate=1405994845498&api=v2

## Postgres docker
docker pull postgres:9.1




