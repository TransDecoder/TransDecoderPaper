# Footprint coverage/support by RNA expression of the source transcript (StringTie TPM).
# usage: awk -f expr_matched_support.awk stringtie.gtf <prefix>.stringtie.orf_support.tsv
# out: class, TPM bin, n, frac with >=10 P-sites, frac frame-supported
BEGIN{FS="\t"}
FNR==NR{ if($3=="transcript"){match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16); match($9,/TPM "[^"]+"/); tpm[t]=substr($9,RSTART+5,RLENGTH-6)+0}; next}
FNR>1{ t=$1; sub(/\.p[0-9]+$/,"",t); x=tpm[t]; b=(x<1?"a<1":x<5?"b1-5":x<20?"c5-20":x<100?"d20-100":"e>=100"); k=$3"\t"b; n[k]++; if($5>=10)c[k]++; s[k]+=$11 }
END{for(k in n) printf "%s\t%d\t%.3f\t%.3f\n",k,n[k],c[k]/n[k],s[k]/n[k]}
