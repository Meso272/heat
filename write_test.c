#include <stdio.h>
int main(int argc, char *argv[]){
    FILE *File = fopen("test.dat", "wb");
    int a=1;
    fseek(File,8,SEEK_SET);
    fwrite(a,sizeof(int),1,File);
    File.close();
    return 0;
}