tfun = function(x = 0){
  if(hasArg(name = 'x')) {
    print('found')
  } else{
    stop('you did not supply x')
  }
}
