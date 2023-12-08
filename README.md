# Hướng Dẫn Sử Dụng SmithWaterman_GPU
## Điều Kiện Set Up Trước Khi Chạy
Để có thể chạy thành công ứng dụng, bạn cần thực hiện các bước set up sau:
Cài đặt CUDA và trình biên dịch nvcc.

## Các Bước Cài Đặt và Chạy Code
### 1. Clone Repository:
- Copy đường link repo.git: https://github.com/tklee-123/SmithWaterman_GPU.git.
- Mở terminal và chạy lệnh sau để sao chép repository về máy local:
```
git clone https://github.com/tklee-123/SmithWaterman_GPU.git
```

### 2. Mở Folder và Di Chuyển tới Thư Mục Web:
- Mở folder vừa clone về trong VSCode.
- Mở terminal trong VSCocde, di chuyển tới thư mục web:
```
cd web
```

### 3. Cài Đặt Dependencies Cho Server:
- Di chuyển vào thư mục SmithWaterman_GPU/web/server và chạy các lệnh sau:
```
cd server
npm install express
node server.js
```
- Sau khi chạy thành công lệnh inline code `node server.js`, terminal sẽ hiển thị:
> Server running on port 8000
- Tiếp theo, mở một terminal mới và di chuyển vào lại thư mục SmithWaterman_GPU/web/server
```
cd web
cd server
npm install react-scripts --save
```

### 4. Quay Trở Lại Thư Mục Web và Cài Đặt Dependencies Cho React App:
```
cd ..
npm install
npm install react-router-dom
npm start
```

### 5. Mở Giao Diện React và Chạy Ứng Dụng:
- Khi giao diện React hiển thị, bạn sẽ thấy Home page (trang chủ) của dự án chứa các thông tin giới thiệu chung. Trên Nav bar, bạn có thể click vào "APPLICATION" hoặc trên hero page click vào nút "Get Start" để di chuyển tới phần ứng dụng.

### 6. Giao Diện Phần Ứng Dụng (APPLICATION):
- Khi giao diện phần ứng dụng hiển thị, bạn sẽ thấy thông tin giới thiệu về ứng dụng và hướng dẫn chạy.

   - Bước 1: Click vào "Choose File", chọn file test có sẵn trong folder của dự án hoặc tự tạo file test theo yêu cầu. File test là file chứa chuỗi ADN với các ký tự là các nucleotit.
  
   - Bước 2: Sau khi hoàn tất chọn file, click vào "Run" để chạy thuật toán. Đợi khoảng 4-5 giây để màn hình hiển thị kết quả.
  
   - Bước 3: Đối chiếu kết quả nhận được với bảng trên màn hình để xác định đặc trưng gen của chuỗi ADN truyền vào.

### 7. Xử Lý Lỗi Hiển Thị:

Nếu có lỗi hiển thị, vui lòng kiểm tra và đảm bảo rằng bạn đã cài đặt thành công môi trường để chạy thuật toán.





