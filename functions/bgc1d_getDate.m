function DateNow = bgc1d_getDate()
% Script to grab the current date

DateNow = datestr(now);
idxreplace=[strfind(DateNow,':'),strfind(DateNow,' '),strfind(DateNow,'-')];
DateNow(idxreplace)='_';
