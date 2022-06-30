#rsync -e 'ssh -p 12898' -avzr root@cloudvm.china-vo.org:/root/lhstat/latest_data ~/PycharmProjects/lhstat/
rsync -avzr ./csst cham@10.3.10.30:~/L1Pipeline/csst/