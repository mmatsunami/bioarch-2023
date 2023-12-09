# 集団サイズ推定ハンズオン

河合　洋介（国立国際医療研究センター）

2023年12月07日作成

## 準備

### サーバへのログイン

まずsshコマンドを用いてサーバにログインしてください。

"user"と"node_name"は割り当てられたものを指定して下さい。

```sh
ssh user@gwa.ddbj.nig.ac.jp
ssh node_name
```

設定ファイルをよみこみます。

```
source /home/bioarchaeology-pg/setup.sh
# Javaが実行できるか確認
java --version
```

### 使用するプログラム

[bcftools](https://github.com/samtools/bcftools)

[PSMC](https://github.com/lh3/psmc)

[hap-IBD](https://github.com/browning-lab/hap-ibd)

[IBDNe](https://faculty.washington.edu/browning/ibdne.html)
今回はスパコンにダウンロード済みのファイルを使います。

#### プログラムの確認
```
bcftools
psmc
java -jar /home/bioarchaeology-pg/kawai/hap-ibd.jar
java -jar /home/bioarchaeology-pg/kawai/ibdne.23Apr20.ae9.jar
```

ヘルプメッセージが表示されることを確認

### 使用するデータ
シミュレーションで作った*Fictus yaponesiae*のゲノムデータ(`FK01.bam`)とフェージング済み日本人ゲノムデータ使います。

#### データの確認
* BAMファイル
```
samtools tview /home/bioarchaeology-pg/kawai/FK01.bam
```
表示がうまくいかない場合は
```
samtools view /home/bioarchaeology-pg/kawai/FK01.bam | less -S
```

* リファレンスファイル
```
less /home/bioarchaeology-pg/kawai/yaponesia_reference.fasta
```

* フェージング済みVCF
```
bcftools view /home/bioarchaeology-pg/kawai/JPT/JPT.phased.chr1.vcf.gz | less -S
```

### 解析ディレクトリの作成
```
# ホームディレクトに移動して作業ディレクトリを作成
cd
mkdir popsize
cd popsize
```

## 実行
### データの準備(PSMC)
* `bcftools`でBAMファイルからバリアント情報を抽出して、コンセンサス配列をfastq形式で抽出する。
```
bcftools mpileup -Ou -f /home/bioarchaeology-pg/kawai/yaponesia_reference.fasta /home/bioarchaeology-pg/kawai/FK01.bam | \
  bcftools call -c - | \
  vcfutils.pl vcf2fq -d 10 -D 100 | \
  gzip -c > FK01.fq.gz
```

>[!NOTE]
>コマンド入力でバックスラッシュ(`\`)に続けてエンターを入力すると改行される。

* 出力ファイルの確認
```
zless FK01.fq.gz
```
>[!NOTE]
>gzip圧縮されたデータは`zless`コマンドで中身を確認できる。


* PSMCに付属するツール(`fq2psmcfa`)で入力形式に変換する
```
fq2psmcfa FK01.fq.gz > FK01.psmcfa
```
* 出力ファイルの確認
```
less FK01.psmcfa
```
### PSMCの実行
```
psmc -t10 -p "10+5*3+4" -o FK01.psmc FK01.psmcfa
```
>[!NOTE]
>-t : 最大の世代数（集団サイズ2N<sub>0</sub>が単位)
>
>-p : 有効集団サイズを推定する時間の分割パターンの指定。"10+5*3+4"は時間を10,3,3,3,3,3,4に分割して、それぞれの区画の有効手段サイズを推定する。
>
>-o : 出力ファイル名（.psmcをつけておく）

* 出力結果の確認
```
less -S FK01.psmc
```

* プロットの出力
```
psmc_plot.pl -u 2.0e-07 -g 1 -x 100 -L -p FK01 FK01.psmc
```
>[!NOTE]
>-u : 突然変異率
>
>-g : 世代時間（年）
>
>-x : 表示する最小の世代
>
>-p : PDFで出力（デフォルトはeps形式)

* 結果の確認

ターミナルで画像ファイルを閲覧することはできないので、ローカルの端末にPDFをダウンロードする。

ローカルのターミナルでscpコマンドで転送する方法
```
scp username-pg@gwa.ddbj.nig.ac.jp:/home/username-pg/psmc/FK01.pdf .
```
WinSCP(windows)やcyberduck(mac)で転送することも可能です。

### IBDNeの実行
#### IBD sharingの計算
* 22番染色体のIBD sharingをhap-ibdで計算します。
```
java -jar hap-ibd.jar gt=/home/bioarchaeology-pg/kawai/JPT/JPT.phased.chr22.vcf.gz \
  map=/home/bioarchaeology-pg/kawai/genetic_map_GRCh38/gmap.chr22.GRCh38.map \
  out=JPT.phased.chr22 nthreads=8
```
>[!NOTE]
>gt= : 入力するフェージング済みのVCFファイル
>
>map= : 遺伝距離が入ったファイル
>
>out= : 出力ファイルのプレフィックス（ibd.txt.gz, hbd.txt.gzが作られる）
>
>nthreads= : 計算に使うCPU数

* 結果の確認
```
zless -S JPT.phased.chr22.ibd.gz
```
|sample1|hap1|sample2|hap2|chr|strat|end|IBD(cM)|
|----|----|----|----|----|----|----|----|
|NA18971|1|NA19057|2|chr22|36125892|36865529|2.198|
|NA18948|2|NA18971|1|chr22|36125892|36947811|2.435|
|NA18943|1|NA19088|2|chr22|36542080|37181105|2.253|
| .. | .. | .. | .. | .. | .. | .. | .. |

IBD長の長い組み合わせを確認する。
```
zcat JPT.phased.chr22.ibd.gz | sort -k8nr | less -S
```

|sample1|hap1|sample2|hap2|chr|strat|end|IBD(cM)|
|----|----|----|----|----|----|----|----|
|NA18945|2|NA19076|1|chr22|10684134|15918351|9.041|
|NA19068|2|NA19080|2|chr22|10684134|15411267|8.981|
|NA18948|1|NA18961|1|chr22|10684134|15285782|8.962|
| .. | .. | .. | .. | .. | .. | .. | .. |

この例ではセントロメア領域に重なっており、偽陽性の可能性が高い。IBDNeは自動的に偽陽性を取り除くがセントロメア領域を除くなど事前処理をした方が良い。

* 1~21番染色体のIBD sharingをまとめて計算する
```
for cn in {1..21};do java -jar /home/bioarchaeology-pg/kawai/hap-ibd.jar gt=/home/bioarchaeology-pg/kawai/JPT/JPT.phased.chr$cn.vcf.gz map=/home/bioarchaeology-pg/kawai/genetic_map_GRCh38/gmap.chr$cn.GRCh38.map nthreads=8 out=JPT.phased.chr$cn; done
```
>[!NOTE]
>この例ではbashのfor構文を使って逐次自動的にコマンドを入力しています。
>遺伝研スパコン一般環境ではジョブスケジューラー(UGE)を使って、一度に複数の計算ノードで実行することが可能です。
>スパコンを使ったゲノム解析ではこのような分散処理を活用することによって解析時間を大幅に減らすことができます。



## 実習
[ヒトゲノムの公共データ](https://sc.ddbj.nig.ac.jp/advanced_guides/advanced_guide_2023#ヒト全ゲノム解析の公共データの再解析データセット)を使ってPSMCで解析する。

遺伝研スパコンの中からは直接アクセス可能です。

`ls /usr/local/shared_data/public-human-genomes/GRCh38/`

* CRAMファイルからPSMCの入力ファイルを作成（時間がないのでハンズオンでは省略）
```
bcftools mpileup -Ou -f Homo_sapiens_assembly38.fasta \
  /usr/local/shared_data/public-human-genomes/GRCh38/1000Genomes/CRAM/NA12878/NA12878.cram | \
  bcftools call -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip -c > CEU_NA12878.fq.gz
fq2psmcfa CEU_NA12878.fq.gz > CEU_NA12878.psmcfa
```
>[!NOTE]
>CRAMは圧縮率の高いマッピングファイルのフォーマットです。
>
>BAMの半分程度の大きさですが、開く時にリファレンス配列を参照する必要があります。

作成済みのpsmcfaファイルが`/home/bioarchaeology-pg/kawai/`に置いてあります。
```
CEU_NA12878.psmcfa
YRI_NA18496.psmcfa
JPT_NA18945.psmcfa
Karitiana_SS6004476.psmcfa
```
* PSMCを実行する

FK01と同じ手順でPSMCを実行してください。

`-p`パラメータで細かく時間を刻むと計算に時間がかかります。

`-p "4+25*2+4+6"`(ヒトゲノムの推奨値)だと２時間ぐらいかかります

