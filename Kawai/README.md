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

### 使用するプログラム

[BEAGLE5.4](http://faculty.washington.edu/browning/beagle/beagle.html)

今回はスパコンにダウンロード済みのファイルを使います。
http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
