# RAD-seq解析ハンズオン

松波　雅俊（琉球大学大学院医学研究科先進ゲノム検査医学講座）

2023年XX月XX日作成

## 準備

### サーバへのログイン

まずsshコマンドを用いてサーバにログインしてください。

"user"と"node_name"は割り当てられたものを指定して下さい。

```sh
ssh user@gwa.ddbj.nig.ac.jp
ssh node_name
```

### 使用するプログラム

```sh
singularity exec -B /lustre8,/home　/usr/local/biotools/s/stacks:2.65--hdcf5f25_0
```

### 使用するデータ

解析には以下のデータを使用します。

シミュレーションによって作成した <i>Fictus yaponesiae</i> 50 個体分（10 個体/集団）の ddRAD-seq データです。

このデータは，EcoRI（切断サイト:G|AATTC）と MseI（切断サイト:T|TAA）という２つの制限酵素でゲノムを断片化し，各ゲノム断片の両端 100bp のペアエンド配列を解読しています。

（*はワイルドカードです）

```sh
#RAD-seq data
/home/bioarchaeology-pg/data/11/rawdata/*.fq.gz

#reference genome
/home/bioarchaeology-pg/data/11/reference/yaponesia_reference.fasta
```


## stacksを使用したチュートリアル

* https://catchenlab.life.illinois.edu/stacks/manual/

まず解析用のディレクトリを作成します。

```sh
mkdir 231210RAD_tutorial
cd 231210RAD_tutorial
```

作業はこのディレクトリで実施します。


### data filtering

最初にstacksのprcocess_radtagsコマンドを用いて、data filteringを実施します。

下記のコマンドでqualtiyの低いリードやアダプター配列を取り除くことができます。

```sh
singularity exec -B /lustre8,/home　/usr/local/biotools/s/stacks:2.65--hdcf5f25_0 \
process_radtags -P -1 /home/bioarchaeology-pg/data/11/rawdata/*_R1.fq.gz -2 /home/bioarchaeology-pg/data/11/rawdata/*_R2.fq.gz -o samples -c -q --renz_1 ecoRI --renz_2 mseI
```
解析がうまくいっていれば、フォルダsamplesに*_R1.1.fq.gz, *_R2.2.fq.gz, *_R1.rem.1.fq.gz, *_R2.rem.2.fq.gzの4つのファイルが作られます。

このうち、*_R1.1.fq.gz, *_R2.2.fq.gzを以降の解析に使用します。

さらに以降のコマンドで動かすためにファイル名を変更する必要があります。

```sh
mv samples/*_R1.1.fq.gz samples/*.1.fq.gz
mv samples/*_R2.2.fq.gz samples/*.2.fq.gz
```
以上がdata filteringで使用するコマンドになります。

しかし、このコマンドを100回も打ち込んで動かすのは大変なので、シェルスクリプトを作成し、slurmでjobを投げることで並列計算させましょう。

まず、シェルスクリプト`process_radtags_all.sh`を作成します。

```sh process_radtags_all.sh
#!/bin/bash
#SBATCH --mem=64G 

SAMPLE=$1

module load singularity/3.8.3

singularity exec -B /lustre8,/home /usr/local/biotools/s/stacks:2.65--hdcf5f25_0 \
process_radtags -P \
-1 /home/bioarchaeology-pg/data/11/rawdata/${SAMPLE}_R1.fq.gz \
-2 /home/bioarchaeology-pg/data/11/rawdata/${SAMPLE}_R2.fq.gz \
-o samples \
-c -q --renz_1 ecoRI --renz_2 mseI
```

このシェルスクリプトでは、コマンドライン上で指定されたサンプル名が${SAMPLE}に代入され、コマンドが実行されます。

slurmではコマンドライン上でループを回すことで並列してjobを投げることができます。

下記の通り、コマンドを打ち込んでください。

```sh
mkdir samples
mkdir process_radtags_all_log

for LOCATION in {FK,OS,SD,SP,TK}; do
  for NUM in {1..20}; do
    if [ ${NUM} -lt 10 ]; then
      ID="${LOCATION}0${NUM}"
    else
      ID="${LOCATION}${NUM}"
    fi
    sbatch -p all -c 1 -n 1 --qos all -o ./process_radtags_all_log/process_radtags_all.${ID}.log \
    -e ./process_radtags_all_log/process_radtags_all.${ID}.error \
    -J ${ID} \
    process_radtags_all.sh ${ID}
  done
done
```
これで100サンプルを同時に処理できます。

標準出力に以下のように表示されれば成功です。

ディレクトリ`samples`にoutput filesができるので、確認して下さい。

ちなみにディレクトリ`process_radtags_all_log`にはjobごとのログファイルが出力されています。

```sh
Submitted batch job XX
Submitted batch job XX
...
```

### de novo assemble

ここではリファレンスゲノムなしでのstacksのパイプラインを紹介します。

#### ustacks

#### sstacks

#### tsv2bam

#### gstacks

#### populations



### with reference genome

ここではリファレンスゲノムありでのstacksのパイプラインを紹介します。

