BEGIN { FS=" " }
NR==FNR { q[$2]=1; next }
$3 >= 700 {
    if ($1 in q) print $1, $2
    if ($2 in q) print $2, $1
}
