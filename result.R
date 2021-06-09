
get_result<-function(norm = 'column', 
					 trans = 'log', 
					 number_cells = 10000,
					 sampleCT='F', 
					 propsample='T', 
					 NormTrans = 'T'){
	
	datasets = c('MacParland2018',
	 'Segerstolpe2016',
	 'Fan2020',
	 'BaronMouse',
	 'Vieira2019',
	 'Darmanis2017',
	 'Shekhar2016',
	 'BaronHuman',
	 'Wilk2020',
	 'Lawlor2017',
	 'Tasic2016',
	 'ZilionisMouse2019',
	 'MaL2019set2',
	 'Zeisel2015',
	 'TravagliniMouseFACS',
	 'Hochgerner1',
	 'Aizarani2019',
	 'DeMicheli2020',
	 'KildisiuteNbPMC',
	 'KildisiuteNbGOSH',
	 'LiH2017_tumor',
	 'Xin2016',
	 'Ximerakis2019',
	 'Hochgerner2',
	 'Elmentaite2020',
	 'Menon2019B',
	 'WangY2019',
	 'Menon2019A',
	 'ZilionisMouse2019fine',
	 'Vieira2019Travaglini',
	 'TravagliniMouseDroplet',
	 'Hrvatin2018',
	 'Wilk2020fine',
	 'Tirosh2016clean',
	 'James2020',
	 'Shekhar2016Menon',
	 'Liao2020',
	 'sims-farber_PBMC',
	 'Vento10X',
	 'KildisiuteAdr',
	 'Tasic2018',
	 'KimN2020')

	bulk_methods = c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elasticNet","lasso","ridge","EPIC",
					 "DSA","ssKL","ssFrobenius","dtangle", "deconf", "proportionsInAdmixture", "EpiDISH","CAMmarker","CDSeq")
    sc_methods = c("MuSiC","BisqueRNA","DWLS","deconvSeq","SCDC","bseqsc","CPM","TIMER")
	all_methods = c(bulk_methods,sc_methods)
	
	ds = c()
	ms = c()
	vs = c()
	for(dataset in datasets){
		for(meth in all_methods){
			if(meth %in% bulk_methods){
				tp = 'bulk'
				name = sprintf('RDS/s.%s.%s.%s.%s.all.%s.%d.none.1.%s.%s.%s.rds', dataset, trans, tp, norm, meth, number_cells, sampleCT, propsample, NormTrans)
			}else{
				tp = 'sc'
				name = sprintf('RDS/s.%s.%s.%s.%s.%s.%s.%d.none.1.%s.%s.%s.rds', dataset, trans, tp, norm, norm, meth, number_cells, sampleCT, propsample, NormTrans)
			}
			if(file.exists(name)){
				x = readRDS(name)
				y = evaluation_metrics(x)$Pearson
			}else{
				y = -1
			}
			ds <- c(ds, dataset)
			ms <- c(ms, meth)
			vs <- c(vs, y)
			
		}
	}
	return(data.frame(cbind(dataset = ds,
      method = ms,
      values = vs)))
	
}

########################
# Plot
colors = c('#1f77b4',
 '#ff7f0e',
 '#279e68',
 '#d62728',
 '#aa40fc',
 '#8c564b',
 '#e377c2',
 '#b5bd61',
 '#17becf',
 '#aec7e8',
 '#ffbb78',
 '#98df8a',
 '#ff9896',
 '#c5b0d5',
 '#c49c94',
 '#f7b6d2',
 '#dbdb8d',
 '#9edae5',
 '#ad494a',
 '#8c6d31',
 '#FFFF00',
 '#1CE6FF',
 '#FF34FF',
 '#FF4A46',
 '#008941',
 '#006FA6',
 '#A30059',
 '#FFDBE5',
 '#7A4900',
 '#0000A6',
 '#63FFAC',
 '#B79762',
 '#004D43',
 '#8FB0FF',
 '#997D87',
 '#5A0007',
 '#809693',
 '#6A3A4C',
 '#1B4400',
 '#4FC601',
 '#3B5DFF',
 '#4A3B53',
 '#FF2F80',
 '#61615A',
 '#BA0900',
 '#6B7900',
 '#00C2A0',
 '#FFAA92',
 '#FF90C9',
 '#B903AA',
 '#D16100',
 '#DDEFFF',
 '#000035',
 '#7B4F4B',
 '#A1C299',
 '#300018',
 '#0AA6D8',
 '#013349',
 '#00846F',
 '#372101',
 '#FFB500',
 '#C2FFED',
 '#A079BF',
 '#CC0744',
 '#C0B9B2',
 '#C2FF99',
 '#001E09',
 '#00489C',
 '#6F0062',
 '#0CBD66',
 '#EEC3FF',
 '#456D75',
 '#B77B68',
 '#7A87A1',
 '#788D66',
 '#885578',
 '#FAD09F',
 '#FF8A9A',
 '#D157A0',
 '#BEC459',
 '#456648',
 '#0086ED',
 '#886F4C',
 '#34362D',
 '#B4A8BD',
 '#00A6AA',
 '#452C2C',
 '#636375',
 '#A3C8C9',
 '#FF913F',
 '#938A81',
 '#575329',
 '#00FECF',
 '#B05B6F',
 '#8CD0FF',
 '#3B9700',
 '#04F757',
 '#C8A1A1',
 '#1E6E00',
 '#7900D7',
 '#A77500',
 '#6367A9',
 '#A05837',
 '#6B002C',
 '#772600',
 '#D790FF',
 '#9B9700',
 '#549E79',
 '#FFF69F',
 '#201625',
 '#72418F',
 '#BC23FF',
 '#99ADC0',
 '#3A2465',
 '#922329',
 '#5B4534',
 '#FDE8DC',
 '#404E55',
 '#0089A3',
 '#CB7E98',
 '#A4E804',
 '#324E72')
                   
plotTopResults<-function(x, n=9, ncols=3){
	library(ggplot2)
	library(DeconRNASeq) 
    parray <- ggplot()
    length(parray) <- n
    i = 1
    for(ts in head(unique(x$tissue), n=n)){
        y = x[x$tissue==ts,]
    parray[[i]]<-ggplot(y, aes(x=observed_values, y=expected_values, color=CT)) + 
        geom_point(alpha=.7)+
        geom_abline(intercept=0, slope=1, colour = "red", size = 1)+ 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
#               text = element_text(size = 20),
              legend.position = "none")+
        scale_color_manual(values=colors)
        i=i+1
    }
    g<-DeconRNASeq::multiplot(plotlist = parray, cols=ncols)
    return(g)
}
plotAllResults<-function(x, n=9, ncols=3, title=''){
	library(ggplot2)
	library(DeconRNASeq) 
    p<-ggplot(x, aes(x=observed_values, y=expected_values, color=CT)) + 
        geom_point(alpha=.7)+
        geom_smooth(method=lm, color='red')+
        ggtitle(title)+ 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              text = element_text(size = 20),
              legend.key = element_rect(colour = "transparent", fill = "white"))+
        scale_color_manual(values=colors)
    return(p)
}

###########################
# Evaluation            
evaluation_metrics<-function(RESULTS, mode=1){
    # return the evaluation metrics
    
    if(mode==1){
        x = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
										   Pearson=cor(observed_values,expected_values) %>% round(.,4))
        y = RESULTS %>% dplyr::group_by(tissue) %>% 
                        dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
                                           Pearson=cor(observed_values,expected_values) %>% round(.,4)) %>% 
                        dplyr::summarise(mRMSE = mean(RMSE), mPearson = mean(Pearson))
        x = cbind(x,y)
    }else{
        x = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values.x-observed_values.y)^2)) %>% round(.,4), 
										   Pearson=cor(observed_values.x,observed_values.y) %>% round(.,4))
        y = RESULTS %>% dplyr::group_by(tissue) %>% 
                        dplyr::summarise(RMSE = sqrt(mean((observed_values.x-expected_values.y)^2)) %>% round(.,4), 
                                           Pearson=cor(observed_values.x,expected_values.y) %>% round(.,4)) %>% 
                        dplyr::summarise(mRMSE = mean(RMSE), mPearson = mean(Pearson))
        x = cbind(x,y)
    }
    return(x)
}
                        
evaluation_table<-function(RESULTS, mode=1){
    # return the evaluation table
    
    if(mode==1){
        x = RESULTS %>% dplyr::group_by(tissue) %>% 
                        dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
                                           Pearson=cor(observed_values,expected_values) %>% round(.,4))
    }else{
        x = RESULTS %>% dplyr::group_by(tissue) %>% 
                        dplyr::summarise(RMSE = sqrt(mean((observed_values.x-expected_values.y)^2)) %>% round(.,4), 
                                           Pearson=cor(observed_values.x,expected_values.y) %>% round(.,4))
    }
    return(x)
}
                        
get_name<-function(params){
    n = ''
    for( i in params){
        j = unlist(strsplit(i, "/"))
        j = gsub('.rds','',dplyr::last(j))
        n = paste0(n, j, sep = ".")
    }
    return(n)
}
