

write_chars_bin(out,"RIFF",4);
write_int(out,0);
write_chars_bin(out,"AVI ",4);

/* kavi->marker=ftell(kavi->out); */
//write_avi_header_chunk(out,&kavi->avi_header,&kavi->stream_header_v,&kavi->stream_format_v);
//write_avi_header_chunk(kavi);
write_chars_bin(out,"LIST",4);
marker1=ftell(out);
write_int(out,0);
write_chars_bin(out,"hdrl",4);
//write_avi_header(out,&kavi->avi_header);

write_chars_bin(out,"avih",4);
marker2=ftell(out);
write_int(out,0);

write_int(out,avi_header->time_delay);
write_int(out,avi_header->data_rate);
write_int(out,avi_header->reserved);
write_int(out,avi_header->flags);
write_int(out,avi_header->number_of_frames);
write_int(out,avi_header->initial_frames);
write_int(out,avi_header->data_streams);
write_int(out,avi_header->buffer_size);
write_int(out,avi_header->width);
write_int(out,avi_header->height);
write_int(out,avi_header->time_scale);
write_int(out,avi_header->playback_data_rate);
write_int(out,avi_header->starting_time);
write_int(out,avi_header->data_length);

t=ftell(out);
fseek(out,marker2,SEEK_SET);
write_int(out,t-marker2-4);
fseek(out,t,SEEK_SET);


write_chars_bin(out,"LIST",4);
sub_marker=ftell(out);
write_int(out,0);
write_chars_bin(out,"strl",4);
write_stream_header(out,&kavi->stream_header_v);
write_stream_format_v(out,&kavi->stream_format_v);

t=ftell(out);
fseek(out,sub_marker,SEEK_SET);
write_int(out,t-sub_marker-4);
fseek(out,t,SEEK_SET);

t=ftell(out);
fseek(out,marker1,SEEK_SET);
write_int(out,t-marker1-4);
fseek(out,t,SEEK_SET);

//write_junk_chunk(out);
int marker3,t;
int r,l,p;
char *junk={ "JUNK IN THE CHUNK! " };

write_chars_bin(out,"JUNK",4);
marker3=ftell(out);
write_int(out,0);

r=4096-ftell(out);
l=strlen(junk);
p=0;

for (t=0; t<r; t++)
{
  putc(junk[p++],out);
  if (p>=l) p=0; 
}

t=ftell(out);
fseek(out,marker3,SEEK_SET);
write_int(out,t-marker3-4);
fseek(out,t,SEEK_SET);





write_chars_bin(out,"LIST",4);
kavi->marker=ftell(out); // size of frames
write_int(out,0);
write_chars_bin(out,"movi",4);
