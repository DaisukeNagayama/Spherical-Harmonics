function getAffine = getAffineFuncHandles()
  affineTranslation = @(tx,ty,tz) [1 0 0 0; 0 1 0 0; 0 0 1 0; tx ty tz 1;];
  affineScaling  = @(sx,sy,sz,cx,cy,cz) affineTranslation(-cx,-cy,-cz)*[sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1;]*affineTranslation(cx,cy,cz);
  affineRotationX = @(rot,cx,cy,cz) affineTranslation(-cx,-cy,-cz)*[1 0 0 0; 0 cos(rot) sin(rot) 0; 0 -sin(rot) cos(rot) 0; 0 0 0 1;]*affineTranslation(cx,cy,cz);
  affineRotationY = @(rot,cx,cy,cz) affineTranslation(-cx,-cy,-cz)*[cos(rot) 0 -sin(rot) 0; 0 1 0 0; sin(rot) 0 cos(rot) 0; 0 0 0 1;]*affineTranslation(cx,cy,cz);
  affineRotationZ = @(rot,cx,cy,cz) affineTranslation(-cx,-cy,-cz)*[cos(rot) sin(rot) 0 0; -sin(rot) cos(rot) 0 0; 0 0 1 0; 0 0 0 1;]*affineTranslation(cx,cy,cz);
  
  getAffine.Translation = affineTranslation;
  getAffine.Scaling = affineScaling;
  getAffine.RotationX = affineRotationX;
  getAffine.RotationY = affineRotationY;
  getAffine.RotationZ = affineRotationZ;
end