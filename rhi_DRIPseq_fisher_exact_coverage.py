import scipy.stats as stats
f1 = open("/Users/ts928/Downloads/rhi_kseek/rhi.DRIPseq.TE2.FET2")
f3 = open("/Users/ts928/Downloads/rhi_kseek/rhi.DRIPseq.TE2.pgreater.FET2.rh",'w')
f1line = f1.readline()
f1line = f1line.rstrip()
f1list = f1line.split('\t')
for n in range(0,13):
        f3.write(f1list[n] + '\t')
f3.write('Ctrl.DRIP.p' + '\t' + 'rhi.DRIP.p' + '\t' +
         'Ctrl.rh.p' + '\t' + 'rhi.rh.p' + '\n')
f3.flush()

for f1_line in f1:
        f1_line = f1_line.rstrip()
        f1_list = f1_line.split('\t')
       # print(f1_list)
        for n in range(0,13):
            f3.write(f1_list[n] + '\t')
        CIP_inbin = int(f1_list[1])
        CIP_notinbin = int(f1_list[2])
        Crh_inbin = int(f1_list[3])
        Crh_notinbin = int(f1_list[4])
        CInp_inbin = int(f1_list[5])
        CInp_notinbin = int(f1_list[6])
        odds_ratio, Cp_val = stats.fisher_exact([[CIP_inbin, CIP_notinbin],
                                                  [CInp_inbin, CInp_notinbin]],
                                                 alternative = 'greater')  
        odds_ratio, Crh_val = stats.fisher_exact([[CIP_inbin, CIP_notinbin],
                                                  [Crh_inbin, Crh_notinbin]],
                                                alternative = 'greater')

        Cp_val = str(Cp_val)
        Crh_val = str(Crh_val)

        rIP_inbin = int(f1_list[7])
        rIP_notinbin = int(f1_list[8])
        rrh_inbin = int(f1_list[9])
        rrh_notinbin = int(f1_list[10])
        rInp_inbin = int(f1_list[11])
        rInp_notinbin = int(f1_list[12])
        odds_ratio, rp_val = stats.fisher_exact([[rIP_inbin, rIP_notinbin],
                                                  [rInp_inbin, rInp_notinbin]],
                                                 alternative = 'greater')
        odds_ratio, rrh_val = stats.fisher_exact([[rIP_inbin, rIP_notinbin],
                                                  [rrh_inbin, rrh_notinbin]],
                                                 alternative = 'greater')
        rp_val = str(rp_val)
        rrh_val = str(rrh_val)
        f3.write(Cp_val + '\t' + rp_val + '\t' +
                 Crh_val + '\t' + rrh_val +'\n')

f1.close()
f3.close()
