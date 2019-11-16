## External libraries

To keep independent:

    find ./external/<library-folder> -type f -exec git update-index --assume-unchanged '{}' \;
