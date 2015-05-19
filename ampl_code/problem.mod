set COLS;
param COL_LB{COLS};
param COL_UB{COLS};
param COL_INT{COLS};
param COL_NAME{COLS} symbolic;

set ROWS;
param ROW_LB{ROWS};
param ROW_UB{ROWS};
param ROW_NAME{ROWS} symbolic;

set ROW_COL dimen 2;
param MAT_VAL{ROW_COL};

set COL_BY_ROW{ r in ROWS} := {(r,c) in ROW_COL};
set ROW_BY_COL{ c in COLS} := {(r,c) in ROW_COL};

set CON_COL := {c in COLS:COL_INT[c]==0};
set INT_COL := COLS diff CON_COL;

set NZ_COL_LB := {c in CON_COL:COL_LB[c]!=0 and COL_LB[c]!=-Infinity};
set NZ_COL_UB := {c in CON_COL:COL_UB[c]!=0 and COL_UB[c]!=+Infinity};

set NZ_ROW_LB := {r in ROWS:ROW_LB[r]!=0 and ROW_LB[r]!=-Infinity};
set NZ_ROW_UB := {r in ROWS:ROW_UB[r]!=0 and ROW_UB[r]!=+Infinity};

param SCALE_COL{COLS} default 0;
param SCALE_ROW{ROWS} default 0;

param SCALE_COL_LB{c in CON_COL} := if c in NZ_COL_LB then COL_LB[c]*(10^SCALE_COL[c]) else COL_LB[c];
param SCALE_COL_UB{c in CON_COL} := if c in NZ_COL_UB then COL_UB[c]*(10^SCALE_COL[c]) else COL_UB[c];

param SCALE_ROW_LB{r in ROWS} := if r in NZ_ROW_LB then ROW_LB[r]*(10^SCALE_ROW[r]) else ROW_LB[r];
param SCALE_ROW_UB{r in ROWS} := if r in NZ_ROW_UB then ROW_UB[r]*(10^SCALE_ROW[r]) else ROW_UB[r];

param LOG_10 := log(10);

var x{c in CON_COL}
, >= SCALE_COL_LB[c]
, <= SCALE_COL_UB[c]
;
var y{c in INT_COL} binary;

minimize obj: x[0];

s.t. ctr{r in ROWS}: 
	SCALE_ROW_LB[r]
	<=
	+sum{c in COL_BY_ROW[r] inter CON_COL} MAT_VAL[r,c]*(10^-SCALE_COL[c])*(10^SCALE_ROW[r])*x[c]
	+sum{c in COL_BY_ROW[r] inter INT_COL} MAT_VAL[r,c]*(10^-SCALE_COL[c])*(10^SCALE_ROW[r])*y[c] 
	<=
	SCALE_ROW_UB[r]
	;

var scale_col{COLS};
var scale_row{ROWS};

minimize scale_criterion : 
	# x- <= x
	+sum{c in NZ_COL_LB}(scale_col[c]+log(abs(COL_LB[c]))/LOG_10)^2
	# x+ >= x
	+sum{c in NZ_COL_UB}(scale_col[c]+log(abs(COL_UB[c]))/LOG_10)^2	
	# b- <= Ax
	+sum{r in NZ_ROW_LB, (r,c) in ROW_COL} (scale_row[r]-scale_col[c]+log(abs(MAT_VAL[r,c]))/LOG_10)^2	
	+sum{r in NZ_ROW_LB}(scale_row[r]+log(abs(ROW_LB[r]))/LOG_10)^2
	# b+ >= Ax
	+sum{r in NZ_ROW_UB}(scale_row[r]+log(abs(ROW_UB[r]))/LOG_10)^2
	+sum{r in NZ_ROW_UB, (r,c) in ROW_COL} (scale_row[r]-scale_col[c]+log(abs(MAT_VAL[r,c]))/LOG_10)^2
	# c-A^t lambda = 0 
	+(scale_col[0])^2
	+sum{(r,c) in ROW_COL}(+log(abs(MAT_VAL[r,c]))/LOG_10+scale_row[r]-scale_col[c])^2
	;
model;
problem lp: x, y, obj, ctr;
problem scaling: scale_col, scale_row, scale_criterion;
