#!/bin/csh

if ( ! -e test.fid ) then
  echo "Error: test.fid is required for this script and the unit test."
  exit 1
endif

nmrPipe -in test.fid | nmrPipe -fn SP -out SP.fid -ov
nmrPipe -in SP.fid | nmrPipe -fn ZF -out SP_ZF.fid -ov
nmrPipe -in SP_ZF.fid | nmrPipe -fn FT -out SP_ZF_FT.ft1 -ov
nmrPipe -in SP_ZF_FT.ft1 | nmrPipe -fn PS -p0 10 -p1 10 -out SP_ZF_FT_PS.ft1 -ov
nmrPipe -in SP_ZF_FT_PS.ft1 | nmrPipe -fn FT -inv -out SP_ZF_FT_PS_FTI.ft1 -ov
