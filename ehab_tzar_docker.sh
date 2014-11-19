docker pull javimarlop/ubuntugis-docker # clone docker image
#docker run --rm -tiv /home/webuser/grassdb:/opt/grassdb javimarlop/ubuntugis-docker # with shared folder
docker run --rm --privileged -ti javimarlop/ubuntugis-docker /bin/bash

useradd -m majavie #--uid xxx
#usermod -a -G fuse majavie
usermod -s /bin/bash majavie

su majavie
cd /home/majavie
mkdir /home/majavie/local0
#mkdir /home/majavie/.ssh
#echo '139.191.147.102' > /home/majavie/.ssh/hosts

sshfs xxx@yyy:/local0/majavie/ /home/majavie/local0

#wget https://tzar-framework.atlassian.net/wiki/download/attachments/4980739/tzar-0.5.4.jar?version=1&modificationDate=1405994845498&api=v2

wget https://tzar-framework.atlassian.net/wiki/download/attachments/4980739/tzar-0.5.5-rc.jar?version=1&modificationDate=1416358965071&api=v2

mv tzar-0.5* tzar.jar

#easy_install -U distribute

#pip install matplotlib

#java -jar tzar.jar execlocalruns http://rdv-framework.googlecode.com/svn/trunk/projects/EMS-paper-eg#lucy_morph

java -jar tzar.jar execlocalruns https://github.com/javimarlop/eHabpy/trunk



###


