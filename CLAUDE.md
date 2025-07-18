# AHC (AtCoder Heuristic Contest) ソリューション

- このソリューションはAHC(AtCoder Heuristic Contest)の問題を1問1プロジェクトで管理しているソリューションです。
- 各プロジェクト内のcppファイルはAtCoderに提出するコードのため、基本的にファイル分割は行えません。

## プロジェクト構造
- 各AHCXXXフォルダ: 1つのコンテスト問題に対応
- in/フォルダ: テスト入力ファイル（0000.txt～）
- out/フォルダ: 出力ファイル
- AHCXXXフォルダ直下のcppファイル: 提出用コード
- *_メモ.md: 各問題の解法メモ

## 開発環境
- Visual Studio プロジェクト（.vcxproj, .sln）
- WSL2環境（Linux）でも動作

## 共通ライブラリ
- Library/AHCTemplate_長期.cpp: 長期コンテスト用テンプレート
- Library/AHCTemplate_短期.cpp: 短期コンテスト用テンプレート

## コーディング規約
- 標準的なインクルード群を使用
- xorshift乱数生成器を使用
- 時間計測機能を組み込み
- using namespace std;を使用

## リファクタリング
- リファクタを行うときは、必ず、デグレが起きないように気をつけ、コンパイルが通るかどうか都度確認しながら、作業をする。
- 変数名、引数名、定数名をリネームするときは、競プロの文脈であることを理解し、短い変数名を使用する。n,m,x,y,u,v,nx,ny,a,b,h,w,tmp,f,flag,cnt,ans,なども積極的に用いる。