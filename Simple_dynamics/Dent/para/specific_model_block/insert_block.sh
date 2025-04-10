for dir in */; do
    folder_name="${dir%/}"
    R_file="${dir}${folder_name}.R"

    # Insert "block = TRUE," after line 335 (which becomes line 336)
    sed -i '' '341a\
block = TRUE,
' "$R_file"

    # Insert "block = TRUE," after line 248 (which becomes line 249)
    sed -i '' '255a\
block = TRUE,
' "$R_file"
done

