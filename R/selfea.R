#' get_statistics_from_dataFrame
#' 
#' @import pwr plyr MASS
#' @description 
#'    A function "get_statistics_from_dataFrame" computes several statistics by reading csv files obtained from input arguments
#' @param
#' df_contrast
#'    A data frame that consists of 'ID' column and expression profile (columns after 'ID' column).
#'    'ID' column should be unique. Column names after 'ID' column should be unique.
#'    Only positive numbers are allowed in expression data. Here is an example.
#' \tabular{rrrrrrrrr}{
#'   \tab ID \tab Y500U100_001 \tab Y500U100_002 \tab Y500U100_003 \tab Y500U100_004 \tab Y500U200_001 \tab Y500U200_002 \tab Y500U200_003 \cr
#'   1 \tab YKL060C \tab 151 \tab 195 \tab 188 \tab 184 \tab 221 \tab 201 \tab 187 \cr
#'   2 \tab YDR155C \tab 154 \tab 244 \tab 237 \tab 232 \tab 190 \tab 187 \tab 215 \cr
#'   3 \tab YOL086C \tab  64 \tab 89 \tab 128 \tab 109 \tab 116 \tab 119 \tab 139 \cr
#'   4 \tab YJR104C \tab 161 \tab 155 \tab 158 \tab 172 \tab 164 \tab 165 \tab 161 \cr
#'   5 \tab YGR192C \tab 157 \tab 161 \tab 173 \tab 175 \tab 177 \tab 164 \tab 176 \cr
#'   6 \tab YLR150W \tab 96 \tab 109 \tab 113 \tab 115 \tab 119 \tab 121 \tab 141 \cr
#'   7 \tab YPL037C \tab 23 \tab 28 \tab 27 \tab 27 \tab 48 \tab 44 \tab 34 \cr
#'   8 \tab YNL007C \tab 53 \tab 58 \tab 64 \tab 63 \tab 66 \tab 70 \tab 78 \cr
#'   9 \tab YBR072W \tab 52 \tab 53 \tab 54 \tab 44 \tab 73 \tab 62 \tab 67 \cr
#'   10 \tab YDR418W_1 \tab 76 \tab 53 \tab 62 \tab 74 \tab 63 \tab 65 \tab 67 \cr
#'   }
#' @param
#' df_group
#'    A data frame that consists of 'Col_Name' and 'Group' columns
#'    This parameter is to match experiment groups to expression profiles of df_contrast.
#'    'Col_Name' should be corresponding to column names of expression profile of df_contrast.
#'    'Group' columns have experiment informaion of columns in expression profile of df_contrast.  Here is an example. See the example of df_contrast together.
#'    \tabular{rrr}{
#'      \tab Col_Name \tab Group \cr
#'    1 \tab Y500U100_001 \tab U100 \cr
#'    2 \tab Y500U100_002 \tab U100 \cr
#'    3 \tab Y500U100_003 \tab U100 \cr
#'    4 \tab Y500U100_004 \tab U100 \cr
#'    5 \tab Y500U200_001 \tab U200 \cr
#'    6 \tab Y500U200_002 \tab U200 \cr
#'    7 \tab Y500U200_003 \tab U200 \cr
#'    }
#' @param
#' padj
#'    Choose one of these c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#'    "fdr" is default option. The option is same to \code{\link[stats]{p.adjust}}.
#' @return
#' A list that consists of the following items:
#'    \tabular{ll}{
#'    $data_table \tab A data frame that consists of ID, GLM Negative Binomial P-value, Cohen's W, GLM Quasi-Poisson P-value, ANOVA with Normal P-value and Cohen's f. \cr
#'    $min_rep \tab Common number of replicates in your group information.  Generally, it is the minimum number of replicates. \cr
#'    $max_rep \tab Maximum number of replicates in your group information. \cr
#'    $nt \tab The number of total experiments in your expression profile. \cr
#'    $ng \tab The number of groups in your group information. \cr
#'    }
#' @examples
#' library(selfea)
#' 
#' ## For this example we will import Gregori data
#' ## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
#' ## An Effect Size Filter Improves the Reproducibility 
#' ## in Spectral Counting-based Comparative Proteomics.
#' ## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')
#' 
#' ## Description:
#' ## Each sample consists in 500ng of standard yeast lisate spiked with 
#' ## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
#' ## The dataset contains a different number of technical replimessagees of each sample
#' 
#' ## import Gregori data
#' data(example_data1)
#' df_contrast <- example_data
#' df_group <- example_group
#' 
#' ## Get statistics through 'get_statistics_from_dataFrame' function
#' list_result <- get_statistics_from_dataFrame(df_contrast,df_group)
#' 
#' ## Get significant features (alpha >= 0.05 and power >= 0.84)
#' significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.84,method='QPF')
#' 
#' @export
#'
get_statistics_from_dataFrame <- function(df_contrast, df_group, padj = 'fdr')
{
  # Check if df_contrast is data frame
  if (!is.data.frame(df_contrast)) {
    stop('df_contrast should be data frame')
  }
  
  # Check if df_group is data frame
  if (!is.data.frame(df_group)) {
    stop('df_group should be data frame')
  }

  # Check if the first column of contrast file is "ID" column
  if (!colnames(df_contrast)[1] == "ID"){
    message("\n'ID' column is not found, so first column is considered as ID column")
    dataset.ID <- df_contrast[,1]
  } else {
    dataset.ID <- df_contrast$ID
  }

  message("Expression data dimension:\t",nrow(df_contrast)," IDs x ",ncol(df_contrast)-1," Columns")
  
  adj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (padj %in% adj_options) {
    message("You chose \"",padj,"\" for p-value adjustment.\n")
  } else {
    stop('Wrong options for p-value adjustment. It should be one of holm, hochberg, hommel, bonferroni, BH, BY, fdr and none')
  }

  # Change factor vector into character vector
  if(is.factor(dataset.ID))
    dataset.ID <- as.character(dataset.ID)
  
  # Check group data frame has two columns
  if (length(colnames(df_group)) != 2) {
    stop("The number of columns of 'Group' from your input should have two columns. Please check the manual.\n")
  } 
  
  # Check if the first column of group file is "Col_Name"
  if (!colnames(df_group)[1] == "Col_Name"){
    message("'Col_Name' column is not found in your Group information, so first column is considered as column names of contrast file")
    colName_group <- df_group[,1]
  } else {
    colName_group <- df_group$Col_Name
  }
  
  # Change factor vector into character vector
  if(is.factor(colName_group))
    colName_group <- as.character(colName_group)
  
  group <- df_group$Group
  # Check if the second column of group file is "Group"
  if (!colnames(df_group)[2] == "Group"){
    message("'Group' column is not found in your Group information, so second column is considered as group names")
    group <- df_group[,2]
  } else {
    group <- df_group$Group
  }
  
  # Change factor vector into character vector
  if(is.factor(group))
    group <- as.character(group)
  
  dataset.expr <- df_contrast[,seq(2,length(df_contrast))]
  colName_expr <- colnames(dataset.expr)
  
  # Check if column names in contrast file and group file are same
  num_col <- length(colName_group)
  for(i in 1:num_col) {
    if (colName_group[i] != colName_expr[i]) {
      stop("There is different column names between contrast file and group file. Check your group file and try again. The order of column names should be same.")
    }
  }
  
  # Count replicates of each groups
  message("Group dimension:")
  exp_groups <- unique(group)
  num_groups <- length(exp_groups)
  num_reps <- c()
  for (j in 1:num_groups){
    group_name <- exp_groups[j]
    num_replicates <- length(group[group==group_name])
    num_reps <- c(num_reps,num_replicates)
    message("\t\t\t\tGroup ",group_name," has ",num_replicates," replicates")
  }

  list_return <- list(data_table = glm_anova(dataset.expr, dataset.ID, group, padj), 
                      min_rep = min(num_reps),
                      max_rep = max(num_reps),
                      nt = num_col,
                      ng = num_groups)
  return (list_return)
}


#' get_statistics_from_file
#' 
#' @import pwr plyr MASS
#' @description 
#'    A function "get_statistics_from_file" computes several statistics by reading csv files obtained from input arguments
#' @param 
#' file_expr
#'    a CSV type file, comma (,) seperated file format, that has unique "ID" at the first column and expression data for the corresponding ID.  Here is an short example.
#'    \tabular{l}{
#'    ID,Y500U100_001,Y500U100_002,Y500U100_003,Y500U100_004,Y500U200_001,Y500U200_002 \cr
#'    YKL060C,151,195,188,184,221,201 \cr
#'    YDR155C,154,244,237,232,190,187 \cr
#'    YOL086C,64,89,128,109,116,119 \cr
#'    }
#' 
#' @param 
#' file_group
#'    a CSV type file, comma (,) seperated file format, that consists of "Col_Name", column names of "file_expr" parameter, and "Group" information of the corresponding column name.
#'    The order of "Col_Name" column have to be same to order of columns in "file_expr".  Here is an example.  See also the example above.
#'    \tabular{l}{
#'    Col_Name,Group \cr
#'    Y500U100_001,U100 \cr
#'    Y500U100_002,U100 \cr
#'    Y500U100_003,U100 \cr
#'    Y500U100_004,U100 \cr
#'    Y500U200_001,U200 \cr
#'    Y500U200_002,U200 \cr
#'    }
#'    
#' @param
#' padj
#'    Choose one of these c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#'    "fdr" is default option.
#' @return 
#' A list that consists of the following items:
#'    \tabular{ll}{
#'    $data_table \tab A data frame that consists of ID, GLM Negative Binomial P-value, Cohen's W, GLM Quasi-Poisson P-value, ANOVA with Normal P-value and Cohen's f. \cr
#'    $min_rep \tab Common number of replicates in your group information.  Generally, it is the minimum number of replicates. \cr
#'    $max_rep \tab Maximum number of replicates in your group information. \cr
#'    $nt \tab The number of total experiments in your expression profile. \cr
#'    $ng \tab The number of groups in your group information. \cr
#'    }
#' @examples
#' library(selfea)
#' 
#' ## For this example we will import Gregori data
#' ## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
#' ## An Effect Size Filter Improves the Reproducibility 
#' ## in Spectral Counting-based Comparative Proteomics.
#' ## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')
#' 
#' ## Description:
#' ## Each sample consists in 500ng of standard yeast lisate spiked with 
#' ## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
#' ## The dataset contains a different number of technical replimessagees of each sample
#' 
#' ## Import Gregori data
#' data(example_data1)
#' df_contrast <- example_data
#' df_group <- example_group
#' 
#' ## Write Gregori data to use 'get_statistics_from_file' function
#' write.csv(df_contrast,"expression.csv",row.names=FALSE)
#' write.csv(df_group,"group.csv",row.names=FALSE)
#' 
#' ## Get statistics
#' list_result <- get_statistics_from_file("expression.csv","group.csv","fdr")
#' 
#' ## Get significant features (alpha >= 0.05 and power >= 0.84)
#' significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.84,method='QPF')
#' 
#' @export
get_statistics_from_file <- function(file_expr = '',file_group = '', padj = 'fdr')
{
  if (file_expr == '') {
    stop('file_expr parameter is empty. Please put valid file path\n')

  } else if (!file.exists(file_expr)){
    stop(file_expr,' does not exist.  Please check the file path\n')
    
  } else if (file_group == '') {
    stop('file_group parameter is empty. Please put valid file path\n')
    
  } else if (!file.exists(file_group)){
    stop(file_group,' does not exist.  Please check the file path\n')
    
  } else {
    message('Expression dataset:\t\t',file_expr)
    message('Group dataset:\t\t\t',file_group)
    
    df_contrast <- read.csv(file_expr, check.names=F)
    df_group <- read.csv(file_group, check.names=F)
  }
  
  ### Use get_statistics_from_dataFrame
  return (get_statistics_from_dataFrame(df_contrast, df_group, padj))
}


#' glm_anova
#' 
#' @import pwr plyr MASS
#' @description 
#'    Calculate P-values from ANOVA using Normal, Quasi-Poisson and Negative Binomial distribution and Cohen's effect sizes
#' @param 
#' dataset.expr
#'    A data frame that has column names for distinguishing experiments and numerical values for expression levels
#' @param
#' dataset.ID
#'    A vector of the obtained expression profile's ID column
#' @param
#' group
#'    A data frame that consists of 'Col_Name' and 'Group' obtained from the user file through get_statistics_from_file.
#' @param
#' padj
#'    Choose one of these c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#'    "fdr" is default option.
#' @return
#' A data frame containing ID, Cohen's W, Cohen's F, Max fold change,
#'    GLM Negative Binomial P-value, GLM Quasi-Poisson P-value and ANOVA with Normal P-value.
glm_anova <- function(dataset.expr, dataset.ID, group, padj ='fdr')
{
  pvals_nor <- c()
  pvals_qp <- c()
  pvals_nb <- c()
  cohen_f <- c()
  cohen_w <- c()
  max_fcs <- c()
  ID_return <- c()
  SC <- NULL
  NR <- nrow(dataset.expr)
  for(i in 1:NR) {
    sc <- as.numeric(dataset.expr[i,])
    if (sum(sc) <= 0) {
      message ("In your data ", dataset.ID[i], " has only zero values. Selfea can't use it, so skip this data.")
    } else {
      ID_return <- c(ID_return, dataset.ID[i])
      df_aov <- data.frame(SC=sc,Run=group)
      
      ### Gaussian
      anv <- aov(SC ~ Run, df_aov)
      sum_anv <- summary(anv)
      pvals_nor <- c(pvals_nor,sum_anv[[1]][["Pr(>F)"]][[1]])
      
      ### ANOVA calculation
      ### http://www.itl.nist.gov/div898/handbook/prc/section4/prc434.htm
      MST <- sum_anv[[1]][["Sum Sq"]][[1]]/sum_anv[[1]][["Df"]][[1]]
      MSE <- sum_anv[[1]][["Sum Sq"]][[2]]/sum_anv[[1]][["Df"]][[2]]
      anv_f <- MST/MSE
      
      ### QuasiPoisson
      model_glm <- glm(SC~Run,data=df_aov,family=quasipoisson(link = "log"))
      model_glm0 <- glm(SC~1,data=df_aov,family=quasipoisson(link = "log"))
      p.value <- anova(model_glm,model_glm0,test="F")$"Pr(>F)"[2]
      pvals_qp <- c(pvals_qp,p.value)
      
      ### GLM using Negative Binomial Dist.
      model_glm <- MASS::glm.nb(SC~Run,data=df_aov)
      model_glm0 <- MASS::glm.nb(SC~1,data=df_aov)
      p.value <- anova(model_glm,model_glm0)$"Pr(Chi)"[2]
      pvals_nb <- c(pvals_nb,p.value)
      
      ### Cohen's f
      numerator <- 0
      N <- nrow(df_aov)
      # http://www.r-tutor.com/elementary-statistics/numerical-measures/variance
      # https://onlinecourses.science.psu.edu/stat501/node/254
      grand_mean <- mean(df_aov$SC)
      df_run2mean <- plyr::ddply(df_aov, "Run", plyr::summarize, mean=mean(SC))
      df_run2nrow <- plyr::ddply(df_aov, "Run", nrow)
      for (j in 1:nrow(df_run2mean)) {
        p = df_run2nrow[j,2]/N
        numerator <- numerator + p*(grand_mean-df_run2mean[j,2])^2
      }
      f <- sqrt(numerator/MSE)
      cohen_f <- c(cohen_f,f)
      
      ### Cohen's w
      num_runs <- length(sc)
      p0 <- c(rep(1/num_runs,num_runs))
      p1 <- df_aov$SC/sum(df_aov$SC)
      cohen_w <- c(cohen_w,pwr::ES.w1(p0,p1))
      
      ### Max fold change
      max_fc <- 0
      df_run2sum <- plyr::ddply(df_aov, "Run", plyr::summarize, sum=sum(SC))
      for (k in 1:nrow(df_run2sum)) {
        for (l in 1:nrow(df_run2sum)) {
          group_numer <- as.character(df_run2sum$Run[k])
          group_denom <- as.character(df_run2sum$Run[l])
          
          denom <- df_run2sum$sum[df_run2sum$Run == group_denom]
          numer <- df_run2sum$sum[df_run2sum$Run == group_numer]
          fc <- numer/denom
          
          # Code for debugging
          #cat(i)
          #cat('\t')
          #cat(df_run2sum$sum[df_run2sum$Run == group_numer])
          #cat('\t')
          #cat(df_run2sum$sum[df_run2sum$Run == group_denom])
          #cat('\t')
          #cat(fc)
          #cat('\n')
          
          if (fc == 'NaN') {
            max_fc <- max_fc
          } else if (fc > max_fc) {
            max_fc <- fc
          }
          
        }
      }
      max_fcs <- c(max_fcs, max_fc)
    }
  }
  
  adj.pval_qp <- p.adjust(pvals_qp,method=padj,n=length(pvals_qp))
  adj.pvals_nb <- p.adjust(pvals_nb,method=padj,n=length(pvals_nb))
  adj.pval_nor <- p.adjust(pvals_nor,method=padj,n=length(pvals_nor))
  
  df_out <- data.frame(Protein=ID_return, Cohens_W=cohen_w, Cohens_F=cohen_f, Max_FC=max_fcs,
                       QP_Pval_adjusted=adj.pval_qp, NB_Pval_adjusted=adj.pvals_nb, Normal_Pval_adjusted=adj.pval_nor)
  return (df_out)
}

#' top_table
#' 
#' @import pwr plyr MASS
#' @description 
#'    Get IDs that pass two filters, p-value and effect-size.  This top_table will make a significant list that is less than p-value and greater than effect-size.  Effect-size are calculated by obtained power level.
#'    This function requires four parameters. ex) top_table(input_data,pvalue=0.05,power_desired=0.84,method='QPF')
#' @param
#' input_list
#'    The list produced by 'get_statistics_from_file' or 'get_statistics_from_dataFrame' function.
#'    See \code{\link{get_statistics_from_file}} and \code{\link{get_statistics_from_dataFrame}} for more information.
#'    It consists of the following items:
#'    \tabular{ll}{
#'        $data_table \tab A data frame that consists of ID, GLM Negative Binomial P-value, Cohen's W, GLM Quasi-Poisson P-value, ANOVA with Normal P-value and Cohen's f. \cr
#'        $min_rep \tab Common number of replicates in your group information.  Generally, it is the minimum number of replicates. \cr
#'        $max_rep \tab Maximum number of replicates in your group information. \cr
#'        $nt \tab The number of total experiments in your expression profile. \cr
#'        $ng \tab The number of groups in your group information. \cr
#'    }
#' @param 
#' pvalue
#'    p-value should be ranged between 0 to 1.
#'    default is 0.05.
#' @param 
#' power_desired
#'    Give the statistical power you desired for output significant list
#' @param 
#' method
#'    Choose statistics method you want to use for making significant list
#'    \tabular{lll}{
#'    1 \tab "QPF" \tab combination of Quasi-Poisson and Cohen's f. Default. \cr
#'    2 \tab "QPFC" \tab combination of Quasi-Poisson and Fold change. \cr
#'    3 \tab "NBW" \tab combination of Negative Binomial and Cohen's w. \cr
#'    4 \tab "NBFC" \tab combination of Negative Binomial and Fold change. \cr
#'    5 \tab "NORF" \tab combination of ANOVA with normal distribution and Cohen's f. \cr
#'    6 \tab "NORFC" \tab combination of ANOVA with normal distribution and Fold change. \cr
#'    }
#' @param
#' FC_threshold
#'    Fold change you want to use. Default is 2.
#' @return 
#' A list containing the follow items.
#'    \tabular{ll}{
#'    top_table \tab a data frame that consists of ID, Cohen's W, Cohen's F, Max fold change, GLM Negative Binomial P-value, GLM Quasi-Poisson P-value and ANOVA with Normal P-value \cr
#'    minimum_cohen_f \tab Minimum Cohen's f found in the top_table \cr
#'    minimum_cohen_w \tab Minimum Cohen's w found in the top_table \cr
#'    minimum_power \tab Minimum statistical power calculated from hypothsis tested for the top_table's entities \cr
#'    alpha \tab Maximum adjusted p-value \cr
#'    num_group \tab The number of groups used for generating the top_table \cr
#'    num_columns \tab The number of columns (samples or experiments) \cr
#'    common_replicates \tab The number of common replicates. \cr
#'    minimum_FC \tab Minimum fold change that was used to generate top_table \cr
#'    }
#' 
#' @examples
#' library(selfea)
#' 
#' ## For this example we will import Gregori data
#' ## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
#' ## An Effect Size Filter Improves the Reproducibility 
#' ## in Spectral Counting-based Comparative Proteomics.
#' ## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')
#' 
#' ## Description:
#' ## Each sample consists in 500ng of standard yeast lisate spiked with 
#' ## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
#' ## The dataset contains a different number of technical replimessagees of each sample
#' 
#' ## import Gregori data
#' data(example_data1)
#' df_contrast <- example_data
#' df_group <- example_group
#' 
#' ## Get statistics through 'get_statistics_from_dataFrame' function
#' list_result <- get_statistics_from_dataFrame(df_contrast,df_group)
#' 
#' ## Get significant features (alpha >= 0.05 and power >= 0.84)
#' significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.84,method='QPF')
#' 
#' @export
#' 
top_table <- function(input_list,pvalue=0.05,power_desired=0.84,method='QPF',FC_threshold = 2)
{
  input_data <- input_list$data_table
  min_rep <- input_list$min_rep
  max_rep <- input_list$max_rep
  nt <- input_list$nt
  ng <- input_list$ng
  
  # check input_data
  if (!is.data.frame(input_data)) {
    stop('top_table function needs data frame produced by get_statistics_from_file function.')
  }
  if (length(grep("Protein",colnames(input_data))) == 0) {
    stop('Protein column is missing.')
  }
  if (length(grep("NB_Pval_adjusted",colnames(input_data))) == 0) {
    stop('NB_Pval_adjusted column is missing.')
  }
  if (length(grep("Cohens_W",colnames(input_data))) == 0) {
    stop('Cohens_W column is missing.')
  }
  if (length(grep("Cohens_F",colnames(input_data))) == 0) {
    stop('Cohens_F column is missing.')
  }
  if (length(grep("Max_FC",colnames(input_data))) == 0) {
    stop('Max_FC column is missing.')
  }
  if (length(grep("QP_Pval_adjusted",colnames(input_data))) == 0) {
    stop('QP_Pval_adjusted column is missing.')
  }
  if (length(grep("Normal_Pval_adjusted",colnames(input_data))) == 0) {
    stop('Normal_Pval_adjusted column is missing.')
  }
  
  # check pvalue
  if (pvalue > 1) {
    stop('P-value threshould cannot be bigger than 1.')
  } else if (pvalue < 0) {
    stop('P-value threshould should be positive number.')
  }
  
  # check power_desired
  if (power_desired > 1) {
    stop('Statistical power cannot be bigger than 1.')
  } else if (power_desired < 0) {
    stop('Statistical power should be positive number.')
  }

  # Make top table by inputted method option
  # "QPF", "QPFC", "NBW", "NBFC", "NORF", "NORFC"
  if (method == "QPF") {
    message ("You chose GLM Quasi-Poisson model and Cohen's f")
    
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, power = power_desired)
    min_f <- as.double(sprintf("%.4f",tmp$f))
    message("Minimum Cohen's f: ",min_f,", when sig.level = ",pvalue," and minimum power = ",power_desired)
    
    df_return <- subset(input_data,input_data$QP_Pval_adjusted < pvalue)
    df_return <- subset(df_return,df_return$Cohens_F >= min_f)
    df_return <- cbind(Protein=as.character(df_return$Protein), 
                       QP_Pval_adjusted=df_return$QP_Pval_adjusted, Cohens_F=df_return$Cohens_F)
    message("The number of significant IDs is ",nrow(df_return))

    # Make a list to return
    list_return <- list(top_table = df_return, 
                      minimum_cohen_f = min_f,
                      minimum_power = power_desired,
                      alpha = pvalue,
                      num_group = ng,
                      common_replicates = min_rep)
    return (list_return)

  } else if (method == "QPFC") {
      message ("You chose GLM Quasi-Poisson model and Fold change")
      message ("Since you chose Fold Change, 'power_desired' will be ignored")
      df_return <- subset(input_data,input_data$QP_Pval_adjusted < pvalue)
      df_return <- subset(df_return, df_return$Max_FC >= FC_threshold)
      
      min_f <- min(df_return$Cohens_F)
      tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, f = min_f)
      message("When sig.level = ",pvalue," with ", FC_threshold, " fold change, ", "minimum power = ",as.double(sprintf("%.4f",tmp$power)))
      message("Min Cohen's F is ",as.double(sprintf("%.4f",min_f)))
      
      df_return <- cbind(Protein=as.character(df_return$Protein), 
                         QP_Pval_adjusted=df_return$QP_Pval_adjusted, Max_FC=df_return$Max_FC)
      message("The number of significant IDs is ",nrow(df_return))

      # Make a list to return
      list_return <- list(top_table = df_return, 
                          minimum_cohen_f = min_f,
                          minimum_power = tmp$power,
                          alpha = pvalue,
                          num_group = ng,
                          minimum_FC = FC_threshold,
                          common_replicates = min_rep)
      return (list_return)
      
  } else if (method == "NBW") {
    message ("You chose GLM Negative Binomial model and Cohen's w")
    tmp <- pwr::pwr.chisq.test(N = nt, df = ng - 1 , sig.level = pvalue, power = power_desired)
    min_w <- as.double(sprintf("%.4f",tmp$w))
    message("Minimum Cohen's w: ",min_w,", when sig.level = ",pvalue," and minimum power = ",power_desired)
    
    df_return <- subset(input_data,input_data$NB_Pval_adjusted < pvalue)
    df_return <- subset(df_return,df_return$Cohens_W >= min_w)
    df_return <- cbind(Protein=as.character(df_return$Protein), 
                       NB_Pval_adjusted=df_return$NB_Pval_adjusted, Cohens_W=df_return$Cohens_W)
    message("The number of significant IDs is ",nrow(df_return))
    
    # Make a list to return
    list_return <- list(top_table = df_return, 
                        minimum_cohen_w = min_w,
                        minimum_power = power_desired,
                        alpha = pvalue,
                        num_group = ng,
                        num_columns = nt)
    return (list_return)
    
  } else if (method == "NBFC") {
    message ("You chose GLM Negative Binomial model and Fold change")
    message ("Since you chose Fold Change, 'power_desired' will be ignored")
    df_return <- subset(input_data,input_data$NB_Pval_adjusted < pvalue)
    df_return <- subset(df_return, df_return$Max_FC >= FC_threshold)
    
    min_w <- min(df_return$Cohens_W)
    tmp <- pwr::pwr.chisq.test(N = nt, df = ng - 1 , sig.level = pvalue, w = min_w)
    message("When sig.level = ",pvalue," with ", FC_threshold, " fold change, ", "minimum power = ",as.double(sprintf("%.4f",tmp$power)))
    message("Min Cohen's W is ",as.double(sprintf("%.4f",min_w)))

    df_return <- cbind(Protein=as.character(df_return$Protein), 
                       NB_Pval_adjusted=df_return$NB_Pval_adjusted, Max_FC=df_return$Max_FC)
    message("The number of significant IDs is ",nrow(df_return))

    # Make a list to return
    list_return <- list(top_table = df_return, 
                        minimum_cohen_w = min_w,
                        minimum_power = tmp$power,
                        alpha = pvalue,
                        num_group = ng,
                        minimum_FC = FC_threshold,
                        num_columns = nt)
    return (list_return)
    
  } else if (method == "NORF") {
    message ("You chose Normal distribution ANOVA and Cohen's f")
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, power = power_desired)
    min_f <- as.double(sprintf("%.4f",tmp$f))
    message("Minimum Cohen's f: ",min_f,", when sig.level = ",pvalue," and minimum power = ",power_desired)
    df_return <- subset(input_data,input_data$Normal_Pval_adjusted < pvalue)
    df_return <- subset(df_return,df_return$Cohens_F >= min_f)
    df_return <- cbind(Protein=as.character(df_return$Protein), 
                       Normal_Pval_adjusted=df_return$Normal_Pval_adjusted, Cohens_F=df_return$Cohens_F)
    message("The number of significant IDs is ",nrow(df_return))

    # Make a list to return
    list_return <- list(top_table = df_return, 
                        minimum_cohen_f = min_f,
                        minimum_power = power_desired,
                        alpha = pvalue,
                        num_group = ng,
                        common_replicates = min_rep)
    return (list_return)    
    
  } else if (method == "NORFC") {
    message ("You chose Normal distribution ANOVA and Fold change")
    message ("Since you chose Fold Change, 'power_desired' will be ignored")
    
    df_return <- subset(input_data,input_data$Normal_Pval_adjusted < pvalue)
    df_return <- subset(df_return, df_return$Max_FC >= FC_threshold)

    min_f <- min(df_return$Cohens_F)
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, f = min_f)
    message("When sig.level = ",pvalue," with ", FC_threshold, " fold change, ", "minimum power = ",as.double(sprintf("%.4f",tmp$power)))
    message("Min Cohen's F is ",as.double(sprintf("%.4f",min_f)))

    df_return <- cbind(Protein=as.character(df_return$Protein), 
                       Normal_Pval_adjusted=df_return$Normal_Pval_adjusted, Max_FC=df_return$Max_FC)
    message("The number of significant IDs is ",nrow(df_return))

    # Make a list to return
    list_return <- list(top_table = df_return, 
                        minimum_cohen_f = min_f,
                        minimum_power = tmp$power,
                        alpha = pvalue,
                        num_group = ng,
                        minimum_FC = FC_threshold,
                        common_replicates = min_rep)
    return (list_return)
    
  } else {
    stop ("Please choose one between 'QPF', 'QPFC', 'NBW', 'NBFC', 'NORF', 'NORFC'.")
  }
}
