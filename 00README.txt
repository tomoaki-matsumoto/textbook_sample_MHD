■ ファイル一覧

  00README.txt		このファイル
* Makefile		コンパイルに使用する。
  TODO.txt		今後のメモ。
  boundary_free.f90	境界条件：自由境界条件
  boundary_periodic.f90	境界条件：周期境界条件
  boundary_windTunnel.F90	境界条件：Wind tunnel with a step 問題
* config.h		設定ファイル
  cure_crash.F90		計算がクラシュしたときの対応
  errornorm.F90		誤差のnoromを計測するコード
  flux_HLLD.F90		HLLD の数値流束
  flux_HLLD_Boris.F90	HLLD (Boris 版)の数値流束。開発途中なので使用不可
  flux_Roe.F90		Roe の数値流束
  flux_RoeM2.F90		RoeM2 の数値流束。開発途中なので使用不可
  flux_scalarAdvection.F90 スカラー線形移流の数値流束
  grid.F90		格子と格子上の物理量を保存する変数を定義。
  init_Orszag_Tang.F90	初期条件：Orszaq Tang 問題
  init_advect.F90		初期条件：移流問題
  init_shocktube.F90	初期条件：衝撃波管問題
  init_wave.F90		初期条件：波の伝播
  init_windTunnel.F90	初期条件：Wind tunnel with a step 問題
  io.F90			入出力
  main.f90		メインプログラム
* parameter.F90		パラメータを設定。次元数、格子点数
* plot1d.gp		Gnuplot： 1次元データの可視化
  plot1d.pro		IDL： 1次元データの可視化
* plot2d.gp		Gnuplot： 2次元データの可視化
  plot2d.pro		IDL： 2次元データの可視化
  readdatap.pro		IDL：データ読み込み
* readplot.py		Python: データの読み込みと可視化
* run.sh			シミュレーション実行スクリプト
  timestep.F90		時間ステップの推進
  util.f90		雑多なルーチン

* 印は課題で編集するファイル
上記以外のファイルは無視してよい。

■ クイックスタート
コンパイルから実行まで

$ ./runs.h

