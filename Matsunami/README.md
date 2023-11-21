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
singularity exec /usr/local/biotools/s/stacks:2.65--hdcf5f25_0
```

### 使用するデータ

解析には以下のデータを使用します。

シミュレーションによって作成した <i>Fictus yaponesiae</i> 50 個体分（10 個体/集団）の ddRAD-seq データです。

このデータは，EcoRI（切断サイト:G|AATTC）と MseI（切断サイト:T|TAA）という２つの制限酵素でゲノムを断片化し，各ゲノム断片の両端 100bp のペアエンド配列を解読しています。

```sh
#RAD-seq data
/home/bioarchaeology-pg/data/11/rawdata/*.fq.gz

#reference genome
/home/bioarchaeology-pg/data/11/reference/yaponesia_reference.fasta
```


## stacksを使用したチュートリアル

* https://catchenlab.life.illinois.edu/stacks/manual/

### de novo assemble


### with reference genome

