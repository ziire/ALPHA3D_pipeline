##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                 TRAITS AND LANDMARKS                 ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## Traits 
traits1 <- c('ISNSL','ISPNS','NSLNA','NABR','NAPNS',      
             'BRLD', 'OPILD','BAOPI')

landmark_pairs1 <- list(c(1,4), c(1,28), c(4,5),c(5,16), c(5,28), 
                        c(16,17),c(20,17), c(21,20))

## only BILATERAL traits
traits <- c('IS_PMdx',    'IS_PMsx',    'PMdx_ZSdx',   'PMsx_ZSsx',   'PMdx_ZIdx',
            'PMsx_ZIsx',  'PMdx_MTdx',  'PMsx_MTsx',   'NSL_ZSdx',    'NSL_ZSsx',
            'NSL_ZIdx',   'NSL_ZIsx',   'BR_PTsx',     'BR_PTdx',     'BR_APETdx',     
            'BR_APETsx',  'PTdx_APETdx','PTsx_APETsx', 'PTdx_BA',     'PTSx_BA',    
            'PTdx_EAMdx', 'PTsx_EAMsx', 'PTdx_ZYGOdx', 'PTsx_ZYGOsx', 'PTdx_TSPdx', 
            'PTsx_TSPsx', 'ZSdx_ZIdx',  'ZSsx_ZIsx',   'ZIdx_MTdx',   'ZIsx_MTsx',   
            'ZIdx_ZYGOdx','ZIsx_ZYGOsx','ZIdx_TSPdx',  'ZIsx_TSPsx',  'MTdx_PNS',    'MTsx_PNS',
            'PNS_APETdx', 'PNS_APETsx', 'APETdx_BA',   'APETsx_BA',   'APETdx_TSdx', 'APETsx_TSsx', 
            'BA_EAMdx',   'BA_EAMsx',   'EAMdx_ZYGOdx','EAMsx_ZYGOsx','ZYGOdx_TSPdx','ZYGOsx_TSPsx',
            'LD_ASdx',    'LD_ASsx',    'PTdx_ASdx',   'PTsx_ASsx',   'JPdx_ASdx',   'JPsx_ASsx')

landmark_pairs <- list(c(1,2),  c(1,3),  c(2,6),  c(3,11), c(2,7),
                       c(3,12), c(2,29), c(3,30), c(4,6),  c(4,11),
                       c(4,7),  c(4,12), c(16,9), c(16,14),c(16,24),
                       c(16,26),c(9,24), c(14,26),c(9,21), c(14,21),
                       c(9,32), c(14,34),c(9,8),  c(14,13),c(9,10),  
                       c(14,15),c(6,7),  c(11,12),c(7,29), c(12,30),
                       c(7,8),  c(12,13),c(7,10), c(12,15),c(29,28),c(30,28),
                       c(28,24),c(28,26),c(24,21),c(26,21),c(24,25),c(26,27),
                       c(21,32),c(21,34),c(32,8), c(34,13),c(8,10), c(13,15),
                       c(17,18),c(17,19),c(9,18), c(14,19),c(22,18),c(23,19))
## for the AVERAGE
name_traits <- c('ISPM_mean',  'PMZS_mean',  'PMZI_mean',  'PMMT_mean',
                 'NSLZS_mean', 'NSLZI_mean', 'BRPT_mean',  'BRAPET_mean', 'PTAPET_mean', 
                 'PTBA_mean',  'PTEAM_mean', 'PTZYGO_mean','PTTSP_mean',  'ZSZI_mean',   
                 'ZIMT_mean',  'ZIZYGO_mean','ZITSP_mean', 'MTPNS_mean',  'PNSAPET_mean', 
                 'APETBA_mean','APETTS_mean','BAEAM_mean', 'EAMZYGO_mean','ZYGOTSP_mean',
                 'LDAS_mean',  'PTAS_mean',  'JPAS_mean')

## new names for traits (with the averaged ones)
order_traits <- c('ISPM','ISPNS','ISNSL','PMZS','PMZI','PMMT','NSLNA','NAPNS','NSLZS',
                  'NSLZI','MTPNS','ZSZI','ZIMT','ZIZYGO','ZITSP','EAMZYGO','ZYGOTSP',
                  'PTZYGO','PTTSP','NABR', 'BRPT','BRAPET','PTAPET','PTBA','PTEAM','LDAS',
                  'BRLD','OPILD','PTAS','PNSAPET','APETBA','APETTS','BAEAM', 'JPAS','BAOPI')

