# 人口推定ハンズオン

河合　洋介（国立国際医療研究センター）

2023年12月06日作成

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

[BEAGLE5.4](http://faculty.washington.edu/browning/beagle/beagle.html)

今回はスパコンにダウンロード済みのファイルを使います。

#### プログラムの確認
`java -jar /home/bioarchaeology-pg/kawai/beagle.22Jul22.46e.jar`

ヘルプメッセージが表示されることを確認

### 使用するデータ
シミュレーションで作った250人の擬似的なヒトゲノムデータ(`yaponeasi_SP_250.vcf.gz`)を使います。

#### データの確認
```
# ファイルの存在の確認
ls -l /home/bioarchaeology-pg/kawai/yaponesia_SP_250.vcf.gz
# ファイルの中身の確認
bcftools view /home/bioarchaeology-pg/kawai/yaponesia_SP_250.vcf.gz | less -S
```

