set -ev

git clone https://github.com/wdecoster/nanotest.git

bash nanotest/make_cram.sh

python nanoget/nanoget.py
