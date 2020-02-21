CENTRAL_SDS<<-3.5
logit<-function(x) {
    x[(x<0)&(x>(-10^((-15))))]=0
    x[(x>1)&((1-x)>(-10^((-15))))]=1
    retval=log(x)
    retval=retval-log(1.0-x)
    retval
}
get_base_trait<-function(trait) {
    gsub("(.+?)(_tech_adj|_gwas_adj|_tech_ntd_adj|_gwas_ntd_adj|_age_sex_meno_adj|_gwas_normalised)?", "\\1",trait)    
}
adjust_variable_for_technical<-function(data_local,data_filters, trait, platform="Coulter", not_time_day=FALSE)
    {

        eval(parse(text=sprintf("data_local$trait=data_local$%s", trait)))

        if(platform=="Coulter")
            {
                naught100_traits=coulter0100_traits
                naught1_traits=coulter01_traits
                positive_traits=coulter_positive_traits

                }
        if(platform=="Sysmex")
            {
                naught100_traits=sysmex0100_traits
                naught1_traits=sysmex01_traits
                positive_traits=sysmex_positive_traits

                }

        if(is.element(get_base_trait(trait), naught100_traits))
            {
                data_local$adjust_scale_trait=logit(data_local$trait/100)
                }else
                    {

                        if(is.element(get_base_trait(trait), naught1_traits))
                            {
                                data_local$adjust_scale_trait=logit(data_local$trait)

                                }else
                                    {
                                        if(is.element(get_base_trait(trait), positive_traits))
                                        {
                                            data_local$adjust_scale_trait=log(data_local$trait)

                                            }
                                        }
                        }

        # throw out data without an instrument
        # note there was a bug in the line below that wja has now corrected 24/3/2016 used to have no effect
                data_local$adjust_scale_trait[is.na(data_local$instrument)]=NA

        # remove infinite values
        data_local$adjust_scale_trait[abs(data_local$adjust_scale_trait)==Inf]=NA
        if(!not_time_day)
            {
                data_local$adjust_scale_trait_drift_adjusted=adjust_variable_for_drift(data_local,data_filters)$adjusted_trait
                }
        else
            {

                data_local$adjust_scale_trait_drift_adjusted=adjust_variable_for_drift_not_time_day(data_local,data_filters)$adjusted_trait

                }
        # identify outlying days
        outlier_day_data=data_local[,c("adjust_scale_trait_drift_adjusted", "day_of_study_acq","instrument")] %>% group_by(day_of_study_acq, instrument) %>% summarise(n=sum(!is.na(adjust_scale_trait_drift_adjusted)),mean=mean(adjust_scale_trait_drift_adjusted,na.rm=TRUE))
        outlier_day_data$z_score=sqrt(outlier_day_data$n)*abs(outlier_day_data$mean-median(data_local$adjust_scale_trait_drift_adjusted,na.rm=TRUE))/mad(data_local$adjust_scale_trait_drift_adjusted,na.rm=TRUE)

        bad_days_data=na.omit(outlier_day_data[(outlier_day_data$n<10)|(outlier_day_data$z_score>8),])
        bad_data=is.element(sprintf("%s:%d",data_local$instrument, data_local$day_of_study_acq),sprintf("%s:%d", bad_days_data$instrument, bad_days_data$day_of_study_acq))
        print(sum(bad_data,na.rm=TRUE))
        data_local$adjust_scale_trait[bad_data]=NA

        if(!not_time_day)
            {
                drift_adjusted_ret_list=adjust_variable_for_drift(data_local,data_filters)
                }
        else
            {
                drift_adjusted_ret_list=adjust_variable_for_drift_not_time_day(data_local,data_filters)
                }

        if(is.element(get_base_trait(trait), naught100_traits))
        {
            drift_adjusted_ret_list$adjusted_trait=100*expit(drift_adjusted_ret_list$adjusted_trait)
            }
        if(is.element(get_base_trait(trait), naught1_traits))
        {
            drift_adjusted_ret_list$adjusted_trait=expit(drift_adjusted_ret_list$adjusted_trait)
            }

        if(is.element(get_base_trait(trait),positive_traits))
        {
            drift_adjusted_ret_list$adjusted_trait=exp(drift_adjusted_ret_list$adjusted_trait)
            }
        drift_adjusted_ret_list
        }

adjust_variable_for_gwas<-function(data_local, trait, study="UK Biobank", pipeline = "", naught100_traits, naught1_traits, positive_traits) {
    trait_file_name = gsub("#", "", gsub("\\/", "", gsub("\\^", "", gsub("%", "", trait))))
    if (file.exists(paste0(Sys.getenv('OUT_DIR'), "/cache/gwas_data_", trait_file_name, "_results_list.Rdata")) & 1 == 0) {
        load(paste0(Sys.getenv('OUT_DIR'), "/cache/gwas_data_", trait_file_name, "_results_list.Rdata"))
        return(return_list)
    }
    
    if(study=="INTERVAL") {
        platform="Sysmex"
    }
    eval(parse(text=sprintf("data_local$trait=as.numeric(data_local$\"%s\")", trait)))

    ##if(platform=="Sysmex") {
      ##naught100_traits=sysmex0100_traits	
      ##naught1_traits=sysmex01_traits
      ##positive_traits=sysmex_positive_traits
    ##}
    print(paste("dc-avfg b", trait, length(na.omit(data_local$trait))))
    if (trait == "DELTA_HGB_tech_adj") {
        data_local$adjust_scale_trait = data_local$trait
    } else if(is.element(get_base_trait(trait), naught100_traits))
    {
        print(paste("dc-avfg b1", trait, "naught100"))
        data_local$adjust_scale_trait=logit(data_local$trait/100)
    } else if(is.element(get_base_trait(trait), naught1_traits))
    {
        print(paste("dc-avfg b1", trait, "naught1"))
        data_local$adjust_scale_trait=logit(data_local$trait)
        
    } else if(is.element(get_base_trait(trait), positive_traits))
    {
        print(paste("dc-avfg b1", trait, "positive"))
        data_local$adjust_scale_trait=log(data_local$trait)
    } else {
        print(paste("TRAIT NOT ASSIGNED TO VECTOR, ARE YOU SURE IT EXISTS IN THE DF?", trait))
        stop(paste("TRAIT NOT ASSIGNED TO VECTOR, ARE YOU SURE IT EXISTS IN THE DF?", trait))
    }
    print(paste("dc-avfg c", trait, length(na.omit(data_local$adjust_scale_trait))))
    
    if(study=="INTERVAL") {
        if (pipeline == "additional") {
            data_local_baseline = data_local[data_local$interval == "baseline",]
            print(paste("dc-avfg c1", trait, length(na.omit(data_local$adjust_scale_trait))))
	    ## check for duplicates in the baseline data (duplicate: the same participant measured twice)
            dup_subs=unique(data_local_baseline$subject_id[duplicated(data_local_baseline$subject_id)])
	    ## loop through each duplicate
            print(paste("dc-avfg c2", trait, length(na.omit(data_local$adjust_scale_trait))))
            for(sub in dup_subs) {	    
	        ## extract all duplicate records
                sub_records=(data_local_baseline$subject_id==sub)&(!is.na(data_local_baseline$adjust_scale_trait))
		## make sure there are atleast 2 duplicate records (otherwise it isn't duplciate)
                print(paste("dc-avfg c3", trait, length(na.omit(data_local$adjust_scale_trait))))
                if(sum(sub_records)<2){next;}
		## keep the earliest measurement
                data_local_baseline$adjust_scale_trait[sub_records&(data_local_baseline$date_time_acq!=min(data_local_baseline$date_time_acq[sub_records],na.rm=TRUE))]=NA
            }
            print(paste("dc-avfg c4", trait, length(na.omit(data_local$adjust_scale_trait))))
            
	    ## check for duplicates in the second measurement data
            data_local_second = data_local[data_local$interval != "baseline",]
            dup_subs=unique(data_local_second$subject_id[duplicated(data_local_second$subject_id)])
            for(sub in dup_subs) {
                sub_records=(data_local_second$subject_id==sub)&(!is.na(data_local_second$adjust_scale_trait))
                if(sum(sub_records)<2){next;}
                data_local_second$adjust_scale_trait[sub_records&(data_local_second$date_time_acq!=min(data_local_second$date_time_acq[sub_records],na.rm=TRUE))]=NA
            }
            
            data_local <- rbind(data_local_baseline, data_local_second)
        } else {
            print(paste("dc-avfg alt1", trait, length(na.omit(data_local$adjust_scale_trait)))) 
            ## if this isn't the additional dataset all records are at baseline, thus the "interval" column doesn't exist
            dup_subs=unique(data_local$subject_id[duplicated(data_local$subject_id)])
            print(paste("dc-avfg alt2", trait, length(na.omit(data_local$adjust_scale_trait)))) 
            for(sub in dup_subs) {
                sub_records=(data_local$subject_id==sub)&(!is.na(data_local$adjust_scale_trait))
                if(sum(sub_records)<2){next;}
                data_local$adjust_scale_trait[sub_records&(data_local$date_time_acq!=min(data_local$date_time_acq[sub_records],na.rm=TRUE))]=NA
            }
        }
        if (pipeline == "additional") {
        prediction_data=data_local[,c("subject_id", "measurement_id", "adjust_scale_trait", "trait","age_acq", "meno", "height", "weight", "pack_years_smoked","smoking_status", "smoking_amount", "drinking_status","alcohol_intake", "interval", "meno_simple")]
        prediction_data$interval = as.factor(prediction_data$interval)
        } else {
        prediction_data=data_local[,c("subject_id", "measurement_id", "adjust_scale_trait", "trait","age_acq", "meno", "height", "weight", "pack_years_smoked","smoking_status", "smoking_amount", "drinking_status","alcohol_intake", "meno_simple")]
        }
                                        #remove subject specific duplicates
    }
    print(paste("dc-avfg c5", trait, length(na.omit(data_local$adjust_scale_trait))))
    prediction_data$central_data=abs(prediction_data$adjust_scale_trait-median(prediction_data$adjust_scale_trait,na.rm=TRUE))<CENTRAL_SDS*mad(prediction_data$adjust_scale_trait,na.rm=TRUE)
    print(paste("dc-avfg c6", trait, length(na.omit(data_local$adjust_scale_trait))))
    prediction_data$height_missing=as.numeric(is.na(prediction_data$height))
    prediction_data$weight_missing=as.numeric(is.na(prediction_data$weight))
    prediction_data$age_acq_missing=as.numeric(is.na(prediction_data$age_acq))
    prediction_data$pack_years_smoked_missing=as.numeric(is.na(prediction_data$pack_years_smoked))
    
    ## this is an arbitrary value since we are putting dummy adjustments in
    prediction_data$height[is.na(prediction_data$height)]=1
    
    prediction_data$log_height=log(as.numeric(prediction_data$height))
    prediction_data$log_weight=log(as.numeric(prediction_data$weight))
    prediction_data$age_acq[is.na(prediction_data$age_acq)]=1
    prediction_data$pack_years_smoked[is.na(prediction_data$pack_years_smoked)]=1
    
    prediction_data$meno[is.na(prediction_data$meno)]="data_missing"
    prediction_data$smoking_status[is.na(prediction_data$smoking_status)]="data_missing"
    prediction_data$smoking_amount[is.na(prediction_data$smoking_amount)]="data_missing"
    prediction_data$drinking_status[is.na(prediction_data$drinking_status)]="data_missing"
    prediction_data$alcohol_intake[is.na(prediction_data$alcohol_intake)]="data_missing"
    prediction_data$meno_collapse=prediction_data$meno
    prediction_data$meno_collapse[prediction_data$meno_collapse=="unsure"|prediction_data$meno_collapse=="hyst"]="data_missing"
    prediction_data$age_acq = as.numeric(prediction_data$age_acq)
    prediction_data$pack_years_smoked = as.numeric(prediction_data$pack_years_smoked)

    prediction_data_nona=na.omit(prediction_data)
    rownames(prediction_data_nona)=prediction_data_nona$measurement_id
    print(paste("dc-avfg GAM start", trait, nrow(prediction_data_nona)))
    if(study=="INTERVAL" & pipeline == "additional") {
        gam_out=gam(adjust_scale_trait~
                        s(age_acq,by=as.factor(meno_simple), k=30, bs="ps")+
                        as.factor(drinking_status)+
                        as.factor(alcohol_intake)+
                        s(pack_years_smoked, k=19, bs="ps")+
                        as.factor(smoking_status)+
                        as.factor(smoking_amount)+
                        as.factor(interval)+
                        age_acq_missing+
                        pack_years_smoked_missing,
                    data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
                             ##s(log_weight, log_height, by=as.factor(meno), k=30, bs="tp")+
                        ##height_missing+
                        ##weight_missing,
    } else if (study == "INTERVAL") {
        gam_out=gam(adjust_scale_trait~
                        s(age_acq,by=as.factor(meno_simple), k=30, bs="ps")+
                        as.factor(drinking_status)+
                        as.factor(alcohol_intake)+
                        s(pack_years_smoked, k=19, bs="ps")+
                        as.factor(smoking_status)+
                        as.factor(smoking_amount)+
                        age_acq_missing+
                        pack_years_smoked_missing,
                    data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
                        ##s(log_weight, log_height, by=as.factor(meno), k=30, bs="tp")+
                        ##height_missing+
                        ##weight_missing,
    } else if (study == "INTERVAL" & pipeline == "additional" & pull_data == TRUE) {
        gam_out=gam(adjust_scale_trait~
                        s(age_acq,by=as.factor(meno_simple), k=30, bs="ps")+
                        as.factor(drinking_status)+
                        as.factor(alcohol_intake)+
                        s(pack_years_smoked, k=19, bs="ps")+
                        as.factor(smoking_status)+
                        as.factor(smoking_amount)+
                        as.factor(interval)+
                        age_acq_missing+
                        s(log_weight, log_height, by=as.factor(meno), k=30, bs="tp")+
                        height_missing+
                        weight_missing+
                        pack_years_smoked_missing,
                    data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
    } else if (study == "INTERVAL" & pull_data == TRUE) {
        gam_out=gam(adjust_scale_trait~
                        s(age_acq,by=as.factor(meno_simple), k=30, bs="ps")+
                        as.factor(drinking_status)+
                        as.factor(alcohol_intake)+
                        s(pack_years_smoked, k=19, bs="ps")+
                        as.factor(smoking_status)+
                        as.factor(smoking_amount)+
                        age_acq_missing+
                        s(log_weight, log_height, by=as.factor(meno), k=30, bs="tp")+
                        height_missing+
                        weight_missing+
                        pack_years_smoked_missing,
                    data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
        print(summary(gam_out)); stop("a")
    }
    print(sprintf("Explained: %f percent of the variance",100*summary(gam_out)$r.sq))
    predicted_vals=predict(gam_out, type="response", newdata=prediction_data_nona, na.action=na.pass)
    
    prediction_data$adjusted_trait=as.vector(
        data_local$adjust_scale_trait-predicted_vals[match(data_local$measurement_id,rownames(predicted_vals))]+
        mean(prediction_data_nona$adjust_scale_trait[prediction_data_nona$central_data], na.rm=TRUE)
    )
    
    ## add line to save sd of data before and after adjustment
    ##df = data.frame(sd_before = sd(data_local$adjust_scale_trait, na.rm=TRUE), sd_after = sd(prediction_data$adjusted_trait, na.rm=TRUE), stringsAsFactors=F)
    ##write.csv(df, paste0(INTERVAL_FBC_DATA_DIR, "/thesis_added/sd_gwas/", make.names(trait), "_gwas_sd.csv"), row.names=F, quote=F)
    ## add line to save gam model summary data
    ##save(list="gam_out", file=paste0(INTERVAL_FBC_DATA_DIR, "/thesis_added/gam_gwas/", make.names(trait), "_gwas_gam.Rdata"))
    
  
    return_list=list()
    return_list$adjusted_trait=prediction_data$adjusted_trait
    return_list$r.sq=summary(gam_out)$r.sq
    return_list$r.sq_full_data=1-var(data_local$adjust_scale_trait-predicted_vals[match(data_local$measurement_id,rownames(predicted_vals))], na.rm=TRUE)/var(data_local$adjust_scale_trait,na.rm=TRUE)
    
    save(list=c("return_list"), file=paste0(Sys.getenv('OUT_DIR'), "/cache/gwas_data_", trait_file_name, "_results_list.Rdata"))
    return(return_list)
}
