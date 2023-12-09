統合生物考古学セミナー
# 集団構造を解析する
木村亮介（琉球大学）

<br><br>
## 0)事前の準備  
WinSCP(Windows)またはCyberduck(Mac)のインストール

WinSCPの遺伝研スパコンでの設定
https://www.genome-sci.jp/old2010-2015/seminar2015/seminar2015_2_Win.pdf

Rのインストール (https://www.r-project.org/)

SplitTree4のインストール (https://software-ab.informatik.uni-tuebingen.de/download/splitstree4/welcome.html)

以下のソフトウェアのインストールはサーバーを使う場合には不要
* Plink (https://www.cog-genomics.org/plink/)
* ADMIXTURE (http://dalexander.github.io/admixture/index.html)
* Plink ( https://www.cog-genomics.org/plink/)
* ADMIXTURE ( http://dalexander.github.io/admixture/index.html)
* EIGENSOFT (https://github.com/DReichLab/EIG.git)
<br><br>
## 1)準備

### 1-1)計算機サーバへの接続
### 1-2)テストファイルディレクトリのコピー
```
cpp -r /home/bioarchaeology-pg/popstruct popstruct
```
### 1-3)ディレクトリの移動
```
cd popstruct
```
### 1-4)ファイルをみる
```
ls
```
### 1-5)ファイルの中身を閲覧
gzファイルの場合zlessを使用。lessで動く場合もあり。
```
zless yaponesia.vcf.gz
``` 
↑↓で移動、qで停止。
yaponesia.vcf.gzを自分のPCにダウンロード
https://drive.google.com/file/d/1VDOa6ynLXHQpAm7w9lRXIioRIhMTaHEM/view?ts=6571a996

<br><br>
## 2)遺伝距離行列を計算して系統ネットワークを作ろう
### 2-1)Rのパッケージをインストール
https://www.bioconductor.org/install/
を参照
SNPRelate, SeqArray, gdsfmt, phangorn
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("SNPRelate")
https://academic.oup.com/bioinformatics/article/28/24/3326/245844?login=true
BiocManager::install("SeqArray")
BiocManager::install("gdsfmt")
install.packages("phangorn")
```
### 2-2)パッケージの呼び出し
```
library("SNPRelate")
library("SeqArray")
library("gdsfmt")
library("phangorn")
```
### 2-2)VCFファイルからGDSフォーマットへの変換
```
vcf.fn <- "yaponesia.vcf.gz"
snpgdsVCF2GDS(vcf.fn,"basic.gds",method="copy.num.of.ref")
snpFromVCFtoGDS <- snpgdsOpen("basic.gds")
```
### 2-3)距離行列の算出
Dij/meanD
```
dissMatrix  <-  snpgdsDiss(snpFromVCFtoGDS, autosome.only=FALSE, remove.monosnp=FALSE, missing.rate=NaN, num.thread=4, verbose=TRUE)
```
### 2-4)NEXUSファイルの作成
```
sampleNumber <- length(read.gdsn(index.gdsn(snpFromVCFtoGDS, "sample.id")))
sampleNameList <- sapply(1:sampleNumber, function(x){strsplit(read.gdsn(index.gdsn(snpFromVCFtoGDS,"sample.id")), split="-M3")[[x]][1]})
outputDist <- dissMatrix$diss
rownames(outputDist) <- sampleNameList
colnames(outputDist) <- sampleNameList
phangorn::writeDist(outputDist, file=paste(vcf.fn, ".nj.nex", sep=""), format="nexus")
```
### 2-5)SplitTree4での図示
yaponesia.vcf.gz.nj.nexをSplitTree4でopen。  
neibour-netネットワーク図が出たことを確認。  
TreesタブからNJを選びApply。 
<br><br>
## 3)ADMIXTURE解析をしよう
### 3-1)aliasの設定
```
alias plink="singularity exec -B /lustre8/home,/home /usr/local/biotools/p/plink\:1.90b6.21--h516909a_0 plink"
alias admixture="singularity exec -B /lustre8/home,/home /usr/local/biotools/a/admixture:1.3.0--0 admixture"
```
### 3-2)VCFファイルをBEDファイルに変換
```
plink --vcf yaponesia.vcf.gz --keep-allele-order --make-bed --out yaponesia
```
bed,bim,famのファイルが生成されたことを確認。
### 3-3)ADMIXTUREのラン
```
admixture yaponesia.bed 2
```
### 3-4)ADMIXTUREの連続ラン(オプションでCVエラーを計算）
``` 
for K in 2 3 4 5;  do admixture --cv yaponesia.bed $K | tee log${K}.out; done  
```

### 3-5）Rで図を描く
yaponesia.4.Qを自分のPCにダウンロード
別の端末を開いて
```
scp XXX-pg@gwa.ddbj.nig.ac.jp:popstruct/yaponesia.4.Q .
```
Rで
```
tbl=read.table("yaponesia.4.Q")
png("k4.png", width = 600, height = 300)
barplot(t(as.matrix(tbl)), col=rainbow(7), xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
```
画像ファイルk4.pngを開く
### 3-6)CVエラーの値
```
grep -h CV log*.out
```
<br><br>
## 4) 主成分分析を実行しよう
### 4-1)aliasの設定
```
alias converf="singularity exec -B /lustre8/home,/home /usr/local/biotools/e/eigensoft\:8.0.0--h2469040_0 convertf"
alias smartpca="singularity exec -B /lustre8/home,/home /usr/local/biotools/e/eigensoft\:8.0.0--h2469040_0 smartpca"
```
### 4-2)BEDファイルをEIGENSTRATフォーマット(GENOファイル)に変換
以下のパラメータファイルを準備

#par.PED.EIGENSTRAT######  
genotypename:    yaponesia.bed  
snpname:         yaponesia.bim  
indivname:       yaponesia.fam  
outputformat:    EIGENSTRAT  
genotypeoutname: yaponesia.geno  
snpoutname:      yaponesia.snp  
indivoutname:    yaponesia.ind  
familynames:     NO  
#########################
```
convertf -p par.PED.EIGENSTRAT
```
### 4-3)集団情報を加えたindファイル（yaponesia.pop.ind）を用意
```
less yaponesia.pop.ind
```  
qで停止。

### 4-4)smartpcaの実行
以下のパラメータファイルを準備

#par.smartpca############  
genotypename:    yaponesia.geno  
snpname:         yaponesia.snp  
indivname:       yaponesia.pop.ind  
evecoutname:     pca_yaponesia.evec  
evaloutname:     pca_yaponesia.eval  
numoutevec:      5  
#########################
```
smartpca -p par.smartpca >pca_yaponesia.log
``` 
### 4-7)Rでプロットの作成
pca_yaponesia.evecを自分のPCにダウンロード
```
fn = "pca_yaponesia.evec"
evec = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "Pop"))
png("PC1vsPC2.png", width = 600, height = 600)
plot(evec$PC1, evec$PC2, col=factor(evec$Pop))
legend("bottomright", legend=levels(factor(evec$Pop)), col=1:length(levels(factor(evec$Pop))), pch=20)
dev.off()
```
<br><br>
### ハンズオンは以上になります。お疲れさまでした。

