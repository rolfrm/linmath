// Based on linmath.h
// https://github.com/datenwolf/linmath.h

// included where needed
// requires math.h

typedef struct{
  union {
    struct{
      float x,y;
    };
    float data[2];
  };
}vec2;


typedef struct{
  union {
    struct{
      float x,y,z;
    };
    float data[3];
    vec2 xy;
  };

}vec3;

typedef struct{
  union {
    struct{
      float x,y,z,w;
    };
    float data[4];
    vec3 xyz;
  };
}vec4;

#define LINMATH_H_OP(n,name, op)			\
  vec##n vec##n##_##name (vec##n a, vec##n const b);

#define LINMATH_H_DEFINE_VEC(n)					\
  LINMATH_H_OP(n,add,+)						\
  LINMATH_H_OP(n,sub,-)						\
  LINMATH_H_OP(n,mul,*)						\
  LINMATH_H_OP(n,div,/)						\
  vec##n vec##n##_scale(vec##n v, float s);			\
  float vec##n##_mul_inner(vec##n a, vec##n b);			\
  float vec##n##_len(vec##n v);					\
  vec##n vec##n##_normalize(vec##n v);				\
  bool vec##n##_compare(vec##n v1, vec##n v2, float eps);	\
  
LINMATH_H_DEFINE_VEC(2)
LINMATH_H_DEFINE_VEC(3)
LINMATH_H_DEFINE_VEC(4)

vec3 vec3_mul_cross(vec3 const a, vec3 const b);
vec3 vec3_reflect(vec3 const v, vec3 const n);
vec4 vec4_mul_cross(vec4 a, vec4 b);
vec4 vec4_reflect(vec4 v, vec4 n);

typedef struct{
  union{
    struct{
      float 
	m00, m01, m02, m03,
	m10, m11, m12, m13,
	m20, m21, m22, m23,
	m30, m31, m32, m33;
    };
    vec4 columns[4];
    float data[4][4];
  };
}mat4;

mat4 mat4_identity();
mat4 mat4_dup(mat4 N);
vec4 mat4_row(mat4 M, int i);
vec4 mat4_col(mat4 M, int i);
mat4 mat4_transpose(mat4 N);
mat4 mat4_add(mat4 a, mat4 b);
mat4 mat4_sub(mat4 a, mat4 b);
mat4 mat4_scale(mat4 a, float k);
mat4 mat4_scale_aniso(mat4 a, float x, float y, float z);
mat4 mat4_mul(mat4 a, mat4 b);
vec4 mat4_mul_vec4( mat4 M, vec4 v);
mat4 mat4_translate(float x, float y, float z);
mat4 mat4_translate_in_place(mat4 M, float x, float y, float z);
mat4 mat4_from_vec3_mul_outer(vec3 a, vec3 b);
mat4 mat4_rotate(mat4 M, float x, float y, float z, float angle);
mat4 mat4_rotate_X(mat4 M, float angle);
mat4 mat4_rotate_Y(mat4 M, float angle);
mat4 mat4_rotate_Z(mat4 M, float angle);
mat4 mat4_invert(mat4 M);
vec3 mat4_orthonormalize(mat4 M);
mat4 mat4_frustum(float l, float r, float b, float t, float n, float f);
mat4 mat4_ortho(float l, float r, float b, float t, float n, float f);
mat4 mat4_perspective(float y_fov, float aspect, float n, float f);
mat4 mat4_look_at(vec3 eye, vec3 center, vec3 up);

typedef vec4 quat;
quat quat_identity();
quat quat_from_axis(vec3 dir, float angle);
quat quat_add(quat a, quat b);
quat quat_sub(quat a, quat b);
quat quat_mul(quat p, quat q);
quat quat_scale(quat v, float s);
float quat_inner_product(quat a, quat b);
quat quat_conj(quat q);
#define quat_normalize vec4_normalize
vec3 quat_mul_vec3(quat q, vec3 v);
mat4 mat4_from_quat(quat q);
mat4 quat_to_mat4(quat q);
mat4 mat4_mul_quat(mat4 M, quat q);
quat quat_from_mat4(mat4 M);
void mat4_print(mat4 mat);
void vec4_print(vec4 v);
void vec3_print(vec3 v);
