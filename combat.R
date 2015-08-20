# read in expression data
dat <- read.csv("expression.csv")

# read in phenotype data
pheno <- read.csv("sample.csv")

# rownames in phenotype should match column names in expression data

# model (e.g. contains 4 covariates)
model = model.matrix(~var1+var2+var3+var4, data = pheno)

# combat adjust the data
combat.expr <- ComBat(dat, pheno$Batch, model)


