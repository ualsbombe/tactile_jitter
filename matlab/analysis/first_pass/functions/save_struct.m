function [] = save_struct(variable, fullpath, overwrite) %#ok<INUSL>
%SAVE_STRUCT Save struct as struct and check whether it should be
%overwritten (automatically checks whether it should be v-7 or v7-3(

two_gigabytes = 2147483648;
temp = whos('variable');
if temp.bytes >= two_gigabytes
    version = '-v7.3';
else
    version = '-v7';
end

if ~exist(fullpath, 'file') || overwrite
    disp(['Saving: ' fullpath]);
    save(fullpath, version, '-struct', 'variable');
else
    disp([fullpath ' already exists; not overwriting']);
end




