#' Gives test statistic value and critical value for one parameter exponential distribution with group size 3
#' @export
#' @param x1 is sample from Exp Population 1.
#' @param x2 is sample from Exp Population 2.
#' @param x3 is sample from Exp Population 3.
#' @param alpha is numeric value (level of significance)
TEST3EXP<-function(x1,x2,x3,alpha)  {



        func1<-function(a,b,c)        #ALRT
                a*log(b/c)

        func2<-function(d,e){       #Heuristic
                d/e
        }

        func3<-function(f,g,h){       #LRT
                (f/h)^g
        }

        k<-3
        nsim<-rep(1,k)
        nsim[1]<-length(x1)
        nsim[2]<-length(x2)
        nsim[3]<-length(x3)
        N<-sum(nsim)


        Value_Heu<-function(x1,x2,x3){
                minimum<-0
                T_tulika<-rep(1,k)
                T_tulika[1]<-mean(x1)-minimum
                T_tulika[2]<-mean(x2)-minimum
                T_tulika[3]<-mean(x3)-minimum
                T_max<-max(T_tulika)
                T_min<-min(T_tulika)
                data.frame(T_max=T_max,T_min=T_min)
        }


        Value_ALRT1<-function(x1,x2,x3){


                minimum<-0

                a1<-nsim[1]
                b1<- N*(mean(x1)-minimum)
                c1<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])
                data.frame(a1=a1,b1=b1,c1=c1)
        }
        Value_ALRT2<-function(x1,x2,x3){
                minimum<-0
                a2<-nsim[2]
                b2<- N*(mean(x2)-minimum)
                c2<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])
                data.frame(a2=a2,b2=b2,c2=c2)
        }
        Value_ALRT3<-function(x1,x2,x3){
                minimum<-0
                a3<-nsim[3]
                b3<- N*(mean(x3)-minimum)
                c3<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])
                data.frame(a3=a3,b3=b3,c3=c3)
        }
        Value_LRT<-function(x1,x2,x3){
                minimum<-0
                T_tulika<-rep(1,k)
                T_tulika[1]<-mean(x1)-minimum
                T_tulika[2]<-mean(x2)-minimum
                T_tulika[3]<-mean(x3)-minimum
                data.frame(alrt1= T_tulika[1],alrt2= T_tulika[2],alrt3= T_tulika[3])
        }
        kruskal_wallis<-function(a,b){
                (12/(N*(N+1)))*(((a^2)/b))-(N+1)
        }


        KW<-function(x1,x2,x3){
                x_bind<-c(x1,x2,x3)
                rank_bind<-rank(x_bind)
                rank_x1<-sum(rank_bind[1:nsim[1]])
                rank_x2<-sum(rank_bind[nsim[1]+1:nsim[2]])
                rank_x3<-sum(rank_bind[nsim[1]+nsim[2]+1:nsim[3]])
                kw_value<-kruskal_wallis(rank_x1,nsim[1])+kruskal_wallis(rank_x2,nsim[2])+kruskal_wallis(rank_x3,nsim[3])
                data.frame(statistic=kw_value)

        }
        KW_Tulika<-KW(x1,x2,x3)$statistic


        Lambda_ALRT_Tulika<- -2*(func1(Value_ALRT1(x1,x2,x3)$a1,Value_ALRT1(x1,x2,x3)$b1,Value_ALRT1(x1,x2,x3)$c1)+func1(Value_ALRT2(x1,x2,x3)$a2,Value_ALRT2(x1,x2,x3)$b2,Value_ALRT2(x1,x2,x3)$c2)+func1(Value_ALRT3(x1,x2,x3)$a3,Value_ALRT3(x1,x2,x3)$b3,Value_ALRT3(x1,x2,x3)$c3))
        Lambda_ALRT_Tulika

        Lambda_LRT_Tulika<- N^N*( func3(Value_LRT(x1,x2,x3)$alrt1, nsim[1],Value_ALRT1(x1,x2,x3)$c1 )* func3(Value_LRT(x1,x2,x3)$alrt2,nsim[2],Value_ALRT2(x1,x2,x3)$c2 )*func3(Value_LRT(x1,x2,x3)$alrt3,nsim[3],Value_ALRT3(x1,x2,x3)$c3))
        Lambda_LRT_Tulika
        Heu_Tulika<- Value_Heu(x1,x2,x3)$T_max/Value_Heu(x1,x2,x3)$T_min

        est_theta= ((nsim[1]*mean(x1))+(nsim[2]*mean(x2))+(nsim[3]*mean(x3)))/N


        boot=1000

        Lambda_ALRT_Tulika_boot<-rep(0,boot)
        Lambda_LRT_Tulika_boot<-rep(0,boot)
        Heu_Tulika_boot<-rep(0,boot)
        KW_Tulika_boot<-rep(0,boot)

        for (j in 1:boot) {
                x1_boot<- rexp(nsim[1],1/est_theta)
                x2_boot<- rexp(nsim[2],1/est_theta)
                x3_boot<- rexp(nsim[3],1/est_theta)
                Lambda_ALRT_Tulika_boot[j]<- -2*(     func1(     Value_ALRT1(x1_boot,x2_boot,x3_boot)$a1,   Value_ALRT1(x1_boot,x2_boot,x3_boot)$b1,    Value_ALRT1(x1_boot,x2_boot,x3_boot)$c1   )+func1(    Value_ALRT2(    x1_boot,x2_boot,x3_boot)$a2,   Value_ALRT2(x1_boot,x2_boot,x3_boot)$b2,     Value_ALRT2(x1_boot,x2_boot,x3_boot)$c2)   +     func1(Value_ALRT3(x1_boot,x2_boot,x3_boot)$a3,Value_ALRT3(x1_boot,x2_boot,x3_boot)$b3,Value_ALRT3(x1_boot,x2_boot,x3_boot)$c3)    )
                Lambda_LRT_Tulika_boot[j]<- N^N*( func3(Value_LRT(x1_boot,x2_boot,x3_boot)$alrt1, nsim[1],Value_ALRT1(x1_boot,x2_boot,x3_boot)$c1 )* func3(Value_LRT(x1_boot,x2_boot,x3_boot)$alrt2,nsim[2],Value_ALRT2(x1_boot,x2_boot,x3_boot)$c2 )*func3(Value_LRT(x1_boot,x2_boot,x3_boot)$alrt3,nsim[3],Value_ALRT3(x1_boot,x2_boot,x3_boot)$c3))
                Heu_Tulika_boot[j]<- Value_Heu(x1_boot,x2_boot,x3_boot)$T_max/Value_Heu(x1_boot,x2_boot,x3_boot)$T_min
                KW_Tulika_boot[j]<- KW(x1_boot,x2_boot,x3_boot)$statistic

        }
        new_lambda_LRT<-Lambda_LRT_Tulika_boot[order(Lambda_LRT_Tulika_boot)]
        C_star_LRT<-new_lambda_LRT[floor(alpha*boot)]
        C_star_LRT
        new_lambda_ALRT<-Lambda_ALRT_Tulika_boot[order(Lambda_ALRT_Tulika_boot)]
        cv_ALRT<-qchisq(p=alpha, df=2, lower.tail=FALSE)
        #C_star_ALRT<-new_lambda_ALRT[floor(0.05*boot)]
        #C_star_ALRT
        new_lambda_Heu<-Heu_Tulika_boot[order(Heu_Tulika_boot)]
        C_star_Heu<-new_lambda_Heu[floor((1-alpha)*boot)]
        new_lambda_KW<-KW_Tulika_boot[order(KW_Tulika_boot)]
        C_star_KW<-new_lambda_KW[floor((1-alpha)*boot)]



        data.frame(Tests=c('LRT','ALRT','U-test','KW' ), test_statistic=c(Lambda_LRT_Tulika,Lambda_ALRT_Tulika, Heu_Tulika, KW_Tulika ), critical_value=c(C_star_LRT, cv_ALRT, C_star_Heu, C_star_KW))
}
