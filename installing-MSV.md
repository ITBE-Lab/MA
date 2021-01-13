# Installing MSV
The below commands installs MSV and its dependencies. As environment, we suggest an Ubuntu 20.04.1 installation.


## Basics

    sudo apt-get -y install build-essential git cmake python3 python3-dev python3-pip zlib1g zlib1g-dev autoconf clang libc++-dev libc++abi-dev
    sudo pip3 install bokeh==1.4.0


## Install PostgreSQL

    sudo sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list'
    wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -
    sudo apt-get update
    sudo apt-get -y install postgresql-12 libpq-dev postgresql-server-dev-12


## Install PostGIS

    sudo apt install postgis postgresql-12-postgis-3
    sudo -u postgres psql
    CREATE EXTENSION postgis;
    \q


## Configure Postgres User

    # log in without password
    sudo -u postgres psql 
    # set password to admin
    \password postgres
in the password prompt enter “admin” as password (without the double quotes)

    # quit
    \q

    # allow md5 connections
    sudo nano /etc/postgresql/12/main/pg_hba.conf

using the editor, replace the line "local all postgres peer" with "local all postgres md5"

    # restart service
    sudo service postgresql restart


## Install MSV as a Python Module

    export PYTHONPATH=$PYTHONPATH:~/buildMA
if you intend to install MSV permanently, you should add the above export statement to your login bash-script (~/.bashrc file). If you omit this addition to your login bash-script, the above export directive gets lost after shell closure.