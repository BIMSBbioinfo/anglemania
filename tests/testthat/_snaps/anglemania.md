# anglemania works correctly with method pearson

    Code
      list_stats(anglemania_object)$mean_zscore[1:10, 1:10]
    Output
                   [,1]        [,2]       [,3]       [,4]        [,5]        [,6]
       [1,]          NA -0.40334364  0.1187222  0.3873909 -0.06951674  0.70424153
       [2,] -0.40334364          NA  0.5872390  0.4761597  0.15648539 -0.70413592
       [3,]  0.11872221  0.58723902         NA  0.1226075  0.78514923 -0.29074216
       [4,]  0.38739086  0.47615972  0.1226075         NA  1.03346348 -0.25130931
       [5,] -0.06951674  0.15648539  0.7851492  1.0334635          NA -0.69883620
       [6,]  0.70424153 -0.70413592 -0.2907422 -0.2513093 -0.69883620          NA
       [7,] -0.07842926  0.32990909  0.1890233 -0.2006030  0.10001995  0.22389845
       [8,] -0.84431736  0.06936535 -0.8740344 -0.3323105 -0.18210264 -1.22160001
       [9,]  0.64960063 -0.76081039  0.4476994  0.4501799  0.18513773  0.66829615
      [10,]  0.45167610  0.41277925  0.4868826 -0.2892201 -0.18434430  0.07768603
                   [,7]         [,8]       [,9]        [,10]
       [1,] -0.07842926 -0.844317362  0.6496006  0.451676105
       [2,]  0.32990909  0.069365354 -0.7608104  0.412779254
       [3,]  0.18902330 -0.874034398  0.4476994  0.486882591
       [4,] -0.20060299 -0.332310502  0.4501799 -0.289220126
       [5,]  0.10001995 -0.182102637  0.1851377 -0.184344302
       [6,]  0.22389845 -1.221600009  0.6682961  0.077686028
       [7,]          NA -0.645410529 -0.7060680  0.106419133
       [8,] -0.64541053           NA  0.1218621  0.009801912
       [9,] -0.70606801  0.121862147         NA -0.143363987
      [10,]  0.10641913  0.009801912 -0.1433640           NA

---

    Code
      str(list_stats(anglemania_object))
    Output
      List of 3
       $ mean_zscore: num [1:9886, 1:9886] NA -0.4033 0.1187 0.3874 -0.0695 ...
       $ sds_zscore : num [1:9886, 1:9886] NA 0.829 0.718 1.804 0.233 ...
       $ sn_zscore  : num [1:9886, 1:9886] NA 0.487 0.165 0.215 0.298 ...

---

    Code
      get_anglemania_genes(anglemania_object)
    Output
        [1] "Gene5069" "Gene7918" "Gene7403" "Gene8577" "Gene6669" "Gene6647"
        [7] "Gene7493" "Gene4746" "Gene5153" "Gene9940" "Gene5236" "Gene2517"
       [13] "Gene7749" "Gene9368" "Gene292"  "Gene3235" "Gene4092" "Gene1206"
       [19] "Gene6233" "Gene1739" "Gene4216" "Gene4346" "Gene8911" "Gene3533"
       [25] "Gene6310" "Gene3628" "Gene4451" "Gene4020" "Gene1418" "Gene5533"
       [31] "Gene5152" "Gene2626" "Gene8040" "Gene1706" "Gene1176" "Gene3827"
       [37] "Gene5751" "Gene4218" "Gene2109" "Gene3717" "Gene5059" "Gene3759"
       [43] "Gene6513" "Gene1760" "Gene8122" "Gene2777" "Gene4938" "Gene5799"
       [49] "Gene8682" "Gene4914" "Gene3692" "Gene7759" "Gene7556" "Gene4814"
       [55] "Gene1239" "Gene2186" "Gene1576" "Gene9135" "Gene7366" "Gene4511"
       [61] "Gene5793" "Gene8564" "Gene2622" "Gene7018" "Gene7076" "Gene7694"
       [67] "Gene4124" "Gene6764" "Gene6956" "Gene959"  "Gene4158" "Gene3994"
       [73] "Gene9300" "Gene8153" "Gene1506" "Gene3809" "Gene5248" "Gene7353"
       [79] "Gene3506" "Gene3487" "Gene9033" "Gene5562" "Gene1328" "Gene7109"
       [85] "Gene2660" "Gene8338" "Gene6400" "Gene1581" "Gene2387" "Gene1632"
       [91] "Gene2782" "Gene2465" "Gene6993" "Gene8324" "Gene4856" "Gene3131"
       [97] "Gene4563" "Gene393"  "Gene1890" "Gene7799" "Gene5216" "Gene346" 
      [103] "Gene2033" "Gene2955" "Gene5677" "Gene223"  "Gene546"  "Gene1914"
      [109] "Gene437"  "Gene1864" "Gene2145" "Gene2473" "Gene8668" "Gene5329"
      [115] "Gene7451" "Gene4520" "Gene4449" "Gene4122" "Gene2364" "Gene3"   
      [121] "Gene4929" "Gene2691" "Gene4831" "Gene3444" "Gene1688" "Gene137" 
      [127] "Gene1502" "Gene2008" "Gene5967" "Gene5325" "Gene7688" "Gene5863"
      [133] "Gene8924" "Gene6510" "Gene1133" "Gene2301" "Gene5512" "Gene7713"
      [139] "Gene6610" "Gene6777" "Gene3378" "Gene9619" "Gene5838" "Gene8892"
      [145] "Gene1361" "Gene6541" "Gene9878" "Gene1921" "Gene5350" "Gene1233"
      [151] "Gene7138" "Gene384"  "Gene1060" "Gene867"  "Gene9463" "Gene122" 
      [157] "Gene4985" "Gene4871" "Gene6731" "Gene6829" "Gene1208" "Gene1320"
      [163] "Gene8050" "Gene9424" "Gene4034" "Gene493"  "Gene2028" "Gene3890"
      [169] "Gene8480" "Gene2920" "Gene9209" "Gene1324" "Gene8975" "Gene9665"
      [175] "Gene2126" "Gene4198" "Gene2906" "Gene2239" "Gene9608" "Gene7194"
      [181] "Gene8984" "Gene4826" "Gene4598" "Gene2236" "Gene3355" "Gene810" 
      [187] "Gene6444" "Gene8727" "Gene1549" "Gene476"  "Gene7330" "Gene9869"
      [193] "Gene9975" "Gene7207" "Gene8584" "Gene2507" "Gene2895" "Gene7216"
      [199] "Gene4762" "Gene8572" "Gene4845" "Gene5445" "Gene3409" "Gene7796"
      [205] "Gene7978" "Gene705"  "Gene2847" "Gene8293" "Gene4196" "Gene4001"
      [211] "Gene2328" "Gene1429" "Gene6938" "Gene4403" "Gene3861" "Gene4593"
      [217] "Gene8056" "Gene4655" "Gene1000" "Gene9227" "Gene6053" "Gene8285"
      [223] "Gene9316" "Gene831"  "Gene1125" "Gene4710" "Gene398"  "Gene3273"
      [229] "Gene9419" "Gene9427" "Gene8295" "Gene5932" "Gene5963" "Gene2069"
      [235] "Gene2963" "Gene5380" "Gene2982" "Gene5185" "Gene8383" "Gene7064"
      [241] "Gene2319" "Gene877"  "Gene7597" "Gene4971" "Gene9074" "Gene2462"
      [247] "Gene2329" "Gene1024" "Gene8938" "Gene567"  "Gene7913" "Gene8573"
      [253] "Gene8894" "Gene1104" "Gene2380" "Gene1456" "Gene1528" "Gene2436"
      [259] "Gene4975" "Gene6147" "Gene2885" "Gene8858" "Gene5280" "Gene8428"
      [265] "Gene9828" "Gene7130" "Gene1852" "Gene7259" "Gene3958" "Gene6413"
      [271] "Gene1946" "Gene5244" "Gene918"  "Gene1773" "Gene5409" "Gene7469"
      [277] "Gene9398" "Gene1805" "Gene5239" "Gene4600" "Gene4693" "Gene5749"
      [283] "Gene695"  "Gene1028" "Gene5342" "Gene3230" "Gene7878" "Gene7503"
      [289] "Gene2189" "Gene6446" "Gene744"  "Gene452"  "Gene7554" "Gene6905"
      [295] "Gene781"  "Gene9023" "Gene3597" "Gene730"  "Gene327"  "Gene4464"
      [301] "Gene1966" "Gene693"  "Gene1314" "Gene1308" "Gene7776" "Gene1916"
      [307] "Gene2908" "Gene4965" "Gene2854" "Gene2803" "Gene8235" "Gene2471"
      [313] "Gene6705" "Gene3354" "Gene4467" "Gene8529" "Gene5307" "Gene5904"
      [319] "Gene8407" "Gene534"  "Gene1673" "Gene1196" "Gene7927" "Gene8514"
      [325] "Gene7028" "Gene9887" "Gene4038" "Gene1079" "Gene840"  "Gene5583"
      [331] "Gene4621" "Gene6361" "Gene3365" "Gene5835" "Gene7402" "Gene1642"
      [337] "Gene8158" "Gene5603" "Gene7813" "Gene162"  "Gene4691" "Gene9166"
      [343] "Gene9391" "Gene4353" "Gene2273" "Gene431"  "Gene4818" "Gene8363"
      [349] "Gene2664" "Gene1486" "Gene2607" "Gene6563" "Gene4884" "Gene4189"
      [355] "Gene4904" "Gene912"  "Gene13"   "Gene4802" "Gene46"   "Gene471" 
      [361] "Gene4340" "Gene5208" "Gene5847" "Gene98"   "Gene8436" "Gene2666"
      [367] "Gene2852" "Gene9938" "Gene4245" "Gene4512" "Gene7034" "Gene9534"
      [373] "Gene1135" "Gene8379" "Gene6855" "Gene7277" "Gene1562" "Gene7772"
      [379] "Gene8290" "Gene8702" "Gene403"  "Gene1010" "Gene3617" "Gene2567"
      [385] "Gene5410" "Gene2105" "Gene1117" "Gene2547" "Gene8006" "Gene1194"
      [391] "Gene391"  "Gene5661" "Gene3535" "Gene1470" "Gene4053" "Gene5639"
      [397] "Gene5126" "Gene5913" "Gene2808" "Gene5979" "Gene6251" "Gene6113"

