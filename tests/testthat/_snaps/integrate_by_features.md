# integrate_by_features integrates Seurat objects correctly using selected features

    Code
      SeuratObject::LayerData(se_integrated)[1:10, 1:10]
    Output
      10 x 10 sparse Matrix of class "dgCMatrix"
    Message
        [[ suppressing 10 column names 'ATGCCAGAACGACT', 'GAACCTGATGAACC', 'TGACTGGATTCTCA' ... ]]
    Output
                                                                         
      GP9     .         .         .         .         . . . . . .        
      CD9     .         .         .         .         . . . . . .        
      TMEM40  .         .         .         .         . . . . . .        
      CA2     .         .         .         .         . . . . . .        
      PF4     .         .         .         .         . . . . . 0.1015229
      ITGA2B  .         .         .         .         . . . . . .        
      GNG11   .         .         .         .         . . . . . .        
      CLU     .         .         .         .         . . . . . .        
      NGFRAP1 0.9263547 0.8172706 1.2047022 0.8493177 . . . . . .        
      TUBB1   0.6901588 0.8740909 0.4129791 0.9259975 . . . . . .        

