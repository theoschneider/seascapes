folder=~/Users/theo/desktop/uni/M/seascapes/

mutselomega --ncat 30 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --until 2000 run_mutsel_bglobin

readmutselomega --every 1 --until 2000 --burnin 1000 --ss run_mutsel_bglobin

