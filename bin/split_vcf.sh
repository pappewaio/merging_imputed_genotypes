awk -vmax="${1}" '
BEGIN {
    FS = "\t";
    file_count = 1;
    max_variants = max;
    max_variants2= max;
    header = "";
    fh=1;
}

# Store header and column header lines
/^#/ { 
  if(header==""){
    header = $0;
  }else{
    header = header "\n" $0; 
  }
  next;
}

# Process data lines
{
    if(fh){
      print header > "split_" file_count ".vcf";
      fh=0;
    }
    if ($2 >= max_variants2) {
        max_variants2=max_variants2+max_variants
        file_count++;
        print header > "split_" file_count ".vcf";
    }
    print $0 >> "split_" file_count ".vcf";
}

'
