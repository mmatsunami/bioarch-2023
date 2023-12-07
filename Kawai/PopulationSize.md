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

今回はスパコンにダウンロード済みのファイルを使います。

#### プログラムの確認
`bcftools`
`psmc`

ヘルプメッセージが表示されることを確認

### 使用するデータ
シミュレーションで作った250人の擬似的なヒトゲノムデータ(`yaponeasi_SP_250.vcf.gz`)を使います。
