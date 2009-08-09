:- module(gff_parser,
          [
           ]).
:- use_foreign_library(foreign(gff_parser4pl)).

%% parse_gff(+File,?GffDoc) is det

%% gff_feature(+GffDoc,?Feature) is nondet
% true if Feature is a feature in GffDoc

%% gff_feature_in_range(+GffDoc,?Feature,+Seq,?Start,?End) is nondet
% true if Feature is a feature in GffDoc overlapping Seq:Start..End


