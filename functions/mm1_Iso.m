  function fmen1 = mm1(var,varIso,Kvar);

   Pvar = max(0,(var+varIso));
   fmen1 = Pvar./(Pvar+Kvar);

  end
